using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
using Rhino;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using System.Linq;
using Rhino.Geometry.Intersect;

namespace VoussoirPlugin03.Components
{
    public class VoussoirCreate : GH_Component
    {
        public VoussoirCreate()
            : base("Springer", "Spr",
                  "Creates a Springer based on the voussoirs closest to the springer line",
                  "Voussoir", "Vault Creation") { }

        public override Guid ComponentGuid => new Guid("EC88F9F2-CD3B-4C41-ADFF-FD189794137C");

        protected override System.Drawing.Bitmap Icon => VoussoirPlugin03.Properties.Resources.springer;

        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddCurveParameter("Springer Line", "SpringerLine", "List of base lines to create vault springers", GH_ParamAccess.list);
            pManager.AddBrepParameter("Voussoirs", "Voussoirs", "Voussoirs to analyse", GH_ParamAccess.tree);
            pManager.AddPlaneParameter("Transversal Planes", "TransversalPlanes", "Planes at each span division", GH_ParamAccess.list);
            pManager.AddNumberParameter("Springer Width", "SpringerWidth", "Distance perpendicular to springer line", GH_ParamAccess.item, 0.3);
        }

        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddBrepParameter("Springers", "S", "Finished Springers", GH_ParamAccess.tree);
            pManager.AddBrepParameter("Voussoirs", "V", "Non transformed voussoirs", GH_ParamAccess.tree);
            //pManager.AddBrepParameter("Log", "L", "All messages generated during execution", GH_ParamAccess.tree);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            //==============================
            // Inputs
            //==============================
            List<Curve> spLines = new List<Curve>();
            if (!DA.GetDataList(0, spLines)) return;
            GH_Structure<GH_Brep> voussoirs = new GH_Structure<GH_Brep>();
            if (!DA.GetDataTree(1, out voussoirs)) return;
            List<GH_Plane> tPlanes = new List<GH_Plane>();
            if (!DA.GetDataList(2, tPlanes)) return;
            double spWidth = 0.3;
            if (!DA.GetData(3, ref spWidth)) return;

            //==============================
            // Remap voussoirs tree: {A;B} -> {B}(A)
            //==============================
            var remappedVoussoirs = new GH_Structure<GH_Brep>();
            foreach (var branchIdx in Enumerable.Range(0, voussoirs.PathCount))
            {
                var oldPath = voussoirs.Paths[branchIdx];
                var branch = voussoirs.Branches[branchIdx];

                if (oldPath.Indices.Length >= 2)
                {
                    int A = oldPath.Indices[0];
                    int B = oldPath.Indices[1];
                    for (int i = 0; i < branch.Count; i++)
                        remappedVoussoirs.Insert(branch[i], new GH_Path(B), A);
                }
                else
                {
                    foreach (var b in branch) remappedVoussoirs.Append(b, oldPath);
                }
            }

            //==============================
            // Extract extrados faces from voussoirs
            //==============================
            var extrados = new GH_Structure<GH_Brep>();
            foreach (var b in remappedVoussoirs.AllData(true).OfType<GH_Brep>())
            {
                if (b.Value.Faces.Count > 2)
                {
                    var extrBrep = b.Value.Faces[1].DuplicateFace(false);
                    extrados.Append(new GH_Brep(extrBrep), new GH_Path(0));
                }
            }

            //==============================
            // Compute closest points & planar distance
            //==============================
            double planarDistance = 0;
            foreach (var sLine in spLines)
            {
                Point3d startXY = new Point3d(sLine.PointAtStart.X, sLine.PointAtStart.Y, 0);
                double minDist = double.MaxValue;

                foreach (var ghBrep in extrados.AllData(true).OfType<GH_Brep>())
                {
                    if (ghBrep.Value.ClosestPoint(sLine.PointAtStart, out Point3d pt, out _, out _, out _, double.MaxValue, out _))
                    {
                        double distXY = new Point3d(pt.X, pt.Y, 0).DistanceTo(startXY);
                        minDist = Math.Min(minDist, distXY);
                    }
                }
                planarDistance = Math.Max(planarDistance, minDist);
            }
            if (spWidth < planarDistance) spWidth = planarDistance;

            //==============================
            // Compute centroid of all vertices
            //==============================
            var allVerts = remappedVoussoirs.AllData(true).OfType<GH_Brep>()
                .Where(b => b.Value != null)
                .SelectMany(b => b.Value.Vertices.Select(v => v.Location))
                .ToList();

            Point3d centroid = new Point3d(
                allVerts.Average(p => p.X),
                allVerts.Average(p => p.Y),
                allVerts.Average(p => p.Z)
            );

            //==============================
            // Initialize output structures
            //==============================
            var almostSpringers = new GH_Structure<GH_Brep>();
            var restVoussoirs = new GH_Structure<GH_Brep>();
            var inVoussoirs = new GH_Structure<GH_Brep>();

            //==============================
            // Main loop over springer lines
            //==============================
            

            foreach (var line in spLines)
            {                              
                // Extrude line to create vertical reference
                var extrusionBrep = Brep.CreateFromSurface(Surface.CreateExtrusion(line, Vector3d.ZAxis * 10));

                // Intersect with voussoirs
                var intersectedVoussoirs = new GH_Structure<GH_Brep>();
                foreach (var branchIdx in Enumerable.Range(0, remappedVoussoirs.Branches.Count))
                {
                    var branch = remappedVoussoirs.Branches[branchIdx];
                    foreach (var ghBrep in branch)
                    {
                        if (ghBrep?.Value == null) continue;
                        if (Rhino.Geometry.Intersect.Intersection.BrepBrep(ghBrep.Value, extrusionBrep, RhinoDoc.ActiveDoc.ModelAbsoluteTolerance, out Curve[] curves, out Point3d[] pts)
                            && (curves.Length > 0 || pts.Length > 0))
                        {
                            intersectedVoussoirs.Append(ghBrep, new GH_Path(branchIdx));
                            inVoussoirs.Append(ghBrep, new GH_Path(branchIdx));
                        }
                    }
                }

                // Compute base surfaces (springerSurfaces)
                var spPlane = new Plane(line.PointAtStart, line.PointAtEnd - line.PointAtStart, Vector3d.ZAxis);
                if (spPlane.DistanceTo(centroid) > 0) spPlane.Flip();

                Line sLine = new Line(line.PointAtStart, line.PointAtEnd);

                var offsetLine = new Line(line.PointAtStart + spPlane.ZAxis * spWidth, line.PointAtEnd + spPlane.ZAxis * spWidth);
                offsetLine.Extend(20, 20);

                var springerLines = new List<Line>();

                foreach (var p in tPlanes)
                {
                    double tA, tB;
                    bool hasA = Intersection.LinePlane(sLine, p.Value, out tA);
                    bool hasB = Intersection.LinePlane(offsetLine, p.Value, out tB);
                    if (hasA && hasB)
                    {
                        springerLines.Add(new Line(sLine.PointAt(tA), offsetLine.PointAt(tB)));
                    }
                }

                var springerSurfacesLocal = new GH_Structure<GH_Brep>();

                for (int i = 0; i < springerLines.Count - 1; i++)
                {
                    var lofts = Brep.CreateFromLoft(new List<Curve> { springerLines[i].ToNurbsCurve(), springerLines[i + 1].ToNurbsCurve() }, Point3d.Unset, Point3d.Unset, LoftType.Normal, false);
                    if (lofts != null && lofts.Length > 0) springerSurfacesLocal.Append(new GH_Brep(lofts[0]), new GH_Path(i));
                }

                // Compute Z lines
                var secondZlineLocal = new GH_Structure<GH_Line>();
                foreach (var branchIdx in Enumerable.Range(0, intersectedVoussoirs.Branches.Count))
                {
                    var branch = intersectedVoussoirs.Branches[branchIdx];
                    var zVals = branch.SelectMany(b => b.Value.Vertices.Select(v => v.Location.Z)).OrderByDescending(z => z).ToList();
                    if (zVals.Count >= 2)
                        secondZlineLocal.Append(new GH_Line(new Line(line.PointAtStart, line.PointAtStart + Vector3d.ZAxis * zVals[1])), new GH_Path(branchIdx));
                }          

                //=========================================
                // 4. Create volumes from springer surfaces and Z-lines
                //=========================================
                GH_Structure<GH_Brep> springerVolumes = new GH_Structure<GH_Brep>();

                for (int s = 0; s < springerSurfacesLocal.Paths.Count; s++)
                {
                    var surf = springerSurfacesLocal.Branches[s][0].Value;
                    var hLine = secondZlineLocal.Branches[s][0].Value;

                    // Compute height vector
                    Vector3d moveVec = Vector3d.ZAxis * hLine.Length;

                    // Duplicate and move surface
                    Brep surf2 = surf.DuplicateBrep();
                    surf2.Translate(moveVec);

                    // Create lofted walls between original and moved surface edges
                    var openBox = surf.Edges
                        .Select((e, i) => Brep.CreateFromLoft(
                            new List<Curve> { e.EdgeCurve, surf2.Edges[i].EdgeCurve },
                            Point3d.Unset, Point3d.Unset, LoftType.Normal, false)?.FirstOrDefault())
                        .Where(b => b != null)
                        .ToList();

                    // Combine all parts: top, bottom, and walls
                    var box = new List<Brep> { surf, surf2 };
                    box.AddRange(openBox);

                    // Join all into a single Brep volume
                    Brep joinedBox = Brep.JoinBreps(box, RhinoDoc.ActiveDoc.ModelAbsoluteTolerance)?.FirstOrDefault();
                    if (joinedBox != null)
                        springerVolumes.Append(new GH_Brep(joinedBox), new GH_Path(s));
                }

                //=========================================
                // 5. Boolean union between volumes and intersected voussoirs
                //=========================================
                for (int i = 0; i < springerVolumes.Paths.Count; i++)
                {
                    var springerBranch = springerVolumes.Branches[i];
                    var voussoirBranch = i < intersectedVoussoirs.Branches.Count ? intersectedVoussoirs.Branches[i] : null;

                    // Collect all Breps to union
                    var toUnion = springerBranch.Select(b => b?.Value)
                                    .Where(b => b != null)
                                    .ToList();

                    if (voussoirBranch != null)
                        toUnion.AddRange(voussoirBranch.Select(b => b?.Value).Where(b => b != null));

                    if (toUnion.Count == 0) continue;

                    try
                    {
                        // Attempt Boolean union
                        var unionResult = Brep.CreateBooleanUnion(toUnion, RhinoDoc.ActiveDoc.ModelAbsoluteTolerance);
                        if (unionResult != null && unionResult.Length > 0)
                        {
                            foreach (var b in unionResult)
                            {
                                b.MergeCoplanarFaces(RhinoDoc.ActiveDoc.ModelAbsoluteTolerance);
                                almostSpringers.Append(new GH_Brep(b), springerVolumes.Paths[i]);
                            }
                        }
                        else
                        {
                            // Fallback: just merge coplanar faces of original Breps
                            foreach (var b in toUnion)
                            {
                                b.MergeCoplanarFaces(RhinoDoc.ActiveDoc.ModelAbsoluteTolerance);
                                almostSpringers.Append(new GH_Brep(b), springerVolumes.Paths[i]);
                            }
                        }
                    }
                    catch (Exception ex)
                    {
                        RhinoApp.WriteLine($"Boolean union failed at path {springerVolumes.Paths[i]}: {ex.Message}");
                        foreach (var b in toUnion)
                        {
                            b.MergeCoplanarFaces(RhinoDoc.ActiveDoc.ModelAbsoluteTolerance);
                            almostSpringers.Append(new GH_Brep(b), springerVolumes.Paths[i]);
                        }
                    }
                }
            }

            for (int i = 0; i < almostSpringers.PathCount; i++)
            {
                
            }

            //==============================/
            // Separate rest voussoirs
            //==============================
            var tol = RhinoDoc.ActiveDoc.ModelAbsoluteTolerance;
            var intersectedData = inVoussoirs.AllData(true).OfType<GH_Brep>()
                .Where(b => b.Value != null)
                .Select(b => (Volume: VolumeMassProperties.Compute(b.Value)?.Volume ?? 0, Centroid: AreaMassProperties.Compute(b.Value)?.Centroid ?? Point3d.Unset))
                .ToList();

            for (int i = 0; i < remappedVoussoirs.Branches.Count; i++)
            {
                var path = remappedVoussoirs.Paths[i];
                foreach (var ghBrep in remappedVoussoirs.Branches[i])
                {
                    if (!(ghBrep is GH_Brep brepGoo) || brepGoo.Value == null) continue;
                    double vol = VolumeMassProperties.Compute(brepGoo.Value)?.Volume ?? 0;
                    Point3d cen = AreaMassProperties.Compute(brepGoo.Value)?.Centroid ?? Point3d.Unset;

                    if (!intersectedData.Any(d => Math.Abs(d.Volume - vol) < tol && d.Centroid.DistanceTo(cen) < tol))
                        restVoussoirs.Append(ghBrep, path);
                }
            }

            //==============================
            // Output
            //==============================
            DA.SetDataTree(0, almostSpringers);
            DA.SetDataTree(1, restVoussoirs);
            //DA.SetDataTree(2, almostSpringers);
        }
    }
}
