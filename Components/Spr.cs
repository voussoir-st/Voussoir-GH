using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
using Rhino;
using Rhino.Display;
using Rhino.Geometry;
using Rhino.Geometry.Intersect;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Security.Cryptography;

namespace VoussoirPlugin03.Components
{
    public static class TreeUtils
    {
        /// <summary>
        /// Partitions each branch of a GH_Structure into sub-branches of a given size.
        /// Equivalent to the Grasshopper "Partition List" component, applied branch-wise.
        /// </summary>
        public static GH_Structure<T> PartitionTree<T>(
            GH_Structure<T> inputTree, 
            int partitionSize)
            where T : IGH_Goo
        {
            if (partitionSize < 1)
                throw new ArgumentException("Partition size must be at least 1.", nameof(partitionSize));
            GH_Structure<T> partitionedTree = new GH_Structure<T>();
            for (int i = 0; i < inputTree.PathCount; i++)
            {
                var path = inputTree.get_Path(i);
                var branch = inputTree.Branches[i];
                for (int j = 0; j < branch.Count; j += partitionSize)
                {
                    int count = Math.Min(partitionSize, branch.Count - j);
                    var chunk = branch.GetRange(j, count);
                    var newPath = path.AppendElement(j / partitionSize);
                    foreach (var item in chunk)
                        partitionedTree.Append(item, newPath);
                }
            }
            return partitionedTree;
        }
        public static GH_Structure<T> DuplicateBranchElements<T>(GH_Structure<T> inputTree)
            where T : class, IGH_Goo, new()
        {
            GH_Structure<T> duplicatedTree = new GH_Structure<T>();

            for (int i = 0; i < inputTree.PathCount; i++)
            {
                var path = inputTree.get_Path(i);
                var branch = inputTree.Branches[i];

                foreach (var item in branch)
                {
                    // Duplicate safely
                    T dup1 = item.Duplicate() as T;
                    T dup2 = item.Duplicate() as T;

                    // Only append non-null items
                    if (dup1 != null) duplicatedTree.Append(dup1, path);
                    if (dup2 != null) duplicatedTree.Append(dup2, path);
                }
            }

            return duplicatedTree;
        }
        /// <summary>
        /// Finds the closest points on curves for a tree of points and a tree of curves.
        /// Mimics the Grasshopper Curve Closest Point component.
        /// </summary>
        public static void CurveClosestPointTree(
            GH_Structure<GH_Point> pointsTree,
            GH_Structure<GH_Curve> curvesTree,
            out GH_Structure<GH_Point> closestPointsTree,
            out GH_Structure<GH_Number> closestParamsTree,
            out GH_Structure<GH_Number> distancesTree)
        {
            closestPointsTree = new GH_Structure<GH_Point>();
            closestParamsTree = new GH_Structure<GH_Number>();
            distancesTree = new GH_Structure<GH_Number>();

            for (int i = 0; i < pointsTree.PathCount; i++)
            {
                var ptPath = pointsTree.get_Path(i);
                var ptBranch = pointsTree.Branches[i];

                // Use the corresponding curve branch
                var curveBranch = curvesTree.Branches[i];

                foreach (var ghPt in ptBranch)
                {
                    Point3d pt;
                    ghPt.CastTo(out pt);

                    double minDist = double.MaxValue;
                    Point3d closestPt = Point3d.Unset;
                    double closestParam = 0.0;

                    foreach (var ghCrv in curveBranch)
                    {
                        Curve crv = ghCrv.Value;
                        
                        if (crv == null) continue;

                        double t;
                        if (crv.ClosestPoint(pt, out t))
                        {
                            Point3d p = crv.PointAt(t);
                            double d = pt.DistanceTo(p);
                            if (d < minDist)
                            {
                                minDist = d;
                                closestPt = p;
                                closestParam = t;
                            }
                        }
                    }

                    closestPointsTree.Append(new GH_Point(closestPt), ptPath);
                    closestParamsTree.Append(new GH_Number(closestParam), ptPath);
                    distancesTree.Append(new GH_Number(minDist), ptPath);
                }
            }
        }
        /// <summary>
        /// Sorts a tree of values according to a tree of numeric keys (branch by branch).
        /// Mimics the Grasshopper Sort List component.
        /// </summary>
        public static GH_Structure<T> SortListTree<T>(
            GH_Structure<GH_Number> keysTree,
            GH_Structure<T> valuesTree) where T : IGH_Goo
        {
            if (keysTree.PathCount != valuesTree.PathCount)
                throw new ArgumentException("Trees must have the same number of branches.");

            GH_Structure<T> sortedTree = new GH_Structure<T>();

            for (int i = 0; i < keysTree.PathCount; i++)
            {
                var path = keysTree.get_Path(i);
                var keys = keysTree.Branches[i].Select(x => x.Value).ToList();
                var values = valuesTree.Branches[i];

                if (keys.Count != values.Count)
                    throw new ArgumentException($"Branch {i} has mismatched lengths.");

                // Create list of tuples
                List<(T value, double key)> pairs = new List<(T, double)>();
                for (int j = 0; j < keys.Count; j++)
                    pairs.Add((values[j], keys[j]));

                // Sort iteratively
                pairs.Sort((a, b) => a.key.CompareTo(b.key));

                foreach (var p in pairs)
                    sortedTree.Append(p.value, path);
            }

            return sortedTree;
        }
        public static GH_Structure<T> TrimTreeDepth<T>(GH_Structure<T> tree) where T : IGH_Goo
        {
            var trimmedTree = new GH_Structure<T>();

            foreach (var path in tree.Paths)
            {
                var branch = tree.get_Branch(path);
                GH_Path newPath;
                Debug.WriteLine($"branch: " + path);
                if (path.Length > 1)
                {
                    int[] parentIndices = path.Indices.Take(path.Length - 1).ToArray();
                    newPath = new GH_Path(parentIndices);
                }
                else
                {
                    newPath = path;
                }
                Debug.WriteLine($"newPath: " + newPath);
                foreach (T item in branch.Cast<T>())
                    trimmedTree.Append(item, newPath);
            }

            return trimmedTree;
        }

    }
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
            pManager.AddCurveParameter("Log", "L", "All messages generated during execution", GH_ParamAccess.tree);
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
            //==============================
            // Build Voussoirs
            //==============================

            GH_Structure<GH_Brep> matchingFacesTree = new GH_Structure<GH_Brep>();

            double tolerance = Rhino.RhinoMath.DefaultAngleTolerance;

            for (int i = 0; i < almostSpringers.PathCount; i++)
            {
                var path = almostSpringers.get_Path(i);
                var branch = almostSpringers.Branches[i];

                foreach (GH_Brep g in almostSpringers.Branches[i])
                {
                    Brep brep = g.Value; // Get the actual Brep
                    
                    // Filter faces whose planes are contained in tPlanes
                    var matchingFaces = brep.Faces
                        .Cast<BrepFace>()
                        .Where(f =>
                        {
                            Plane facePlane;
                            if (f.TryGetPlane(out facePlane))
                            {
                                // check against all planes in tPlanes
                                foreach (var ghPlane in tPlanes) // tPlanes is List<GH_Plane>
                                {
                                    Plane p = ghPlane.Value; // Extract the Rhino.Geometry.Plane from GH_Plane
                                    // Normals are parallel
                                    if (facePlane.Normal.IsParallelTo(p.Normal, tolerance) != 0)
                                    {
                                        // Distance from origin of one plane to the other plane is small
                                        if (Math.Abs(facePlane.DistanceTo(p.Origin)) < Rhino.RhinoMath.ZeroTolerance)
                                        {
                                            return true;
                                        }
                                    }
                                }
                            }
                            return false;
                        })
                        .Select(f => new GH_Brep(f.DuplicateFace(false)))
                        .ToList();

                    foreach (var face in matchingFaces)
                    {
                        matchingFacesTree.Append(face, path);
                    }
                }
            }

            int partitionSize = 2;
            GH_Structure<GH_Brep> partitionedfacesTree = TreeUtils.PartitionTree(matchingFacesTree, partitionSize);
            partitionedfacesTree.Graft(GH_GraftMode.GraftAll);

            GH_Structure<GH_Curve> edgesTree = new GH_Structure<GH_Curve>();
            
            Debug.WriteLine($"Original tree has {edgesTree.PathCount} paths");
            

            GH_Structure <GH_Point> edgesPointsTree = new GH_Structure<GH_Point>();

            for (int i = 0; i < partitionedfacesTree.PathCount; i++)
            {
                var path = partitionedfacesTree.get_Path(i);
                var branch = partitionedfacesTree.Branches[i];

                foreach (var ghSurface in branch)
                {
                    Brep brep = ghSurface.Value;
                    Curve[] edgeCurves = brep.DuplicateEdgeCurves(false);

                    var joinedCurves = Curve.JoinCurves(edgeCurves);
                    foreach (var curve in joinedCurves)
                    {
                        edgesTree.Append(new GH_Curve(curve), path);
                    }
                    edgesPointsTree.AppendRange(
                        brep.Vertices.Select(v => new GH_Point(v.Location)),
                        path
                    );
                }
            }

            GH_Structure<GH_Curve> trimmedEdgesTree = TreeUtils.TrimTreeDepth(edgesTree);
            Debug.WriteLine($"Trimmed tree has {trimmedEdgesTree.PathCount} paths");

            var dupSpLines = new GH_Structure<GH_Curve>();

            for (int i = 0; i < matchingFacesTree.PathCount; i++)
            {
                dupSpLines.AppendRange(spLines.Select(c => new GH_Curve(c)));
            }
            
            GH_Structure<GH_Curve> dupSpLinesPTree = TreeUtils.PartitionTree(dupSpLines, 2);
            dupSpLinesPTree.Simplify(GH_SimplificationMode.CollapseAllOverlaps);
            dupSpLinesPTree.Graft(GH_GraftMode.GraftAll);

            GH_Structure<GH_Curve> dupSpLinesPTree2 = TreeUtils.DuplicateBranchElements(dupSpLinesPTree);
            dupSpLinesPTree2.Graft(GH_GraftMode.GraftAll);

            GH_Structure<GH_Point> closestPointsTree = new GH_Structure<GH_Point>();
            GH_Structure<GH_Number> closestParamsTree = new GH_Structure<GH_Number>();
            GH_Structure<GH_Number> distancesTree = new GH_Structure<GH_Number>();

            TreeUtils.CurveClosestPointTree(edgesPointsTree, dupSpLinesPTree2, out closestPointsTree, out closestParamsTree, out distancesTree);
            GH_Structure<GH_Point> sortedPointsTree = TreeUtils.SortListTree(distancesTree, closestPointsTree);
            GH_Structure<GH_Point> firstPointsTree = new GH_Structure<GH_Point>();

            for (int i = 0; i < sortedPointsTree.PathCount; i++)
            {
                var path = sortedPointsTree.get_Path(i);
                var branch = sortedPointsTree.Branches[i];

                if (branch.Count > 0)
                    firstPointsTree.Append(branch[0], path);
            }
            
            GH_Structure<GH_Point> trimmed = TreeUtils.TrimTreeDepth(firstPointsTree);

            GH_Structure<GH_Point> a = new GH_Structure<GH_Point>();
            GH_Structure<GH_Number> pr = new GH_Structure<GH_Number>();
            GH_Structure<GH_Number> bc = new GH_Structure<GH_Number>();

            TreeUtils.CurveClosestPointTree(trimmed, trimmedEdgesTree, out a, out pr, out bc);


            //==============================
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
            DA.SetDataTree(1, inVoussoirs);
            DA.SetDataTree(2, trimmedEdgesTree);
        }
    }  
}
