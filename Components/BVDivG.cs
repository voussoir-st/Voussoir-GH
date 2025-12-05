using Components; // Ensure this is present to access Vault utilities (PolylineUtils, Utils, etc.)
using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
using Rhino;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using System.Linq;
using VoussoirPlugin03.Components;
using VoussoirPlugin03.Properties;

namespace Components
{
    public class IrregularVaultDivisionComponent : GH_Component
    {
        public IrregularVaultDivisionComponent()
          : base(
                "Barrel Vault Division - Grid",
                "BVDivG",
                "Divides a vault defined by two arcs into spanwise and lengthwise voussoirs. Accepts a tree of surfaces (each branch one or more surfaces).",
                "Voussoir",
                "Core Geometry"
                )
        { }

        public override Guid ComponentGuid => new Guid("BC98D9F2-CD3B-4C41-ADFF-FD189794437C");

        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                return VoussoirPlugin03.Properties.Resources.VDiv;
            }
        }

        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            // Now accepts a tree of surfaces
            pManager.AddSurfaceParameter("Vault Surface", "VaultSurface", "Surface(s) that define the vault. Provide as a tree where each branch corresponds to one vault (1 or multiple surfaces).", GH_ParamAccess.tree);
            pManager.AddIntegerParameter(
                "Transversal Divisions", "TransversalDivisions",
                "Number of voussoirs in the vault's span\nMinimum: 3",
                GH_ParamAccess.tree, 12);

            pManager.AddIntegerParameter(
                "Longitudinal Divisions", "LongitudinalDivisions",
                "Number of voussoirs in the vault's length\nMinimum: 1",
                GH_ParamAccess.tree, 8);

        }

        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddPlaneParameter("Intrados Planes", "IP", "Intrados planar vault panels (tree matching input branches).", GH_ParamAccess.tree);
            pManager.AddPlaneParameter("Division planes", "DP", "Planes of each Voussoir Contact Surface (tree matching input branches).", GH_ParamAccess.tree);
            pManager.HideParameter(1);
            pManager.AddPlaneParameter("Transversal planes", "TP", "Planes at each span division (tree matching input branches).", GH_ParamAccess.list);
            pManager.HideParameter(2);
            pManager.AddPlaneParameter("Longitudinal planes", "LP", "Planes at each length division (tree matching input branches).", GH_ParamAccess.list);
            pManager.HideParameter(3);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // Input tree of surfaces
            var vaultTree = new GH_Structure<GH_Surface>();
            if (!DA.GetDataTree(0, out vaultTree)) return;

            // Read span and length divisions as trees
            GH_Structure<GH_Integer> spanTree;
            GH_Structure<GH_Integer> lengthTree;

            DA.GetDataTree(1, out spanTree);
            DA.GetDataTree(2, out lengthTree);

            // Flatten span fallback list
            List<int> spanFallback = new List<int>();
            foreach (var p in spanTree.Paths)
            {
                foreach (var ghInt in spanTree[p])
                    spanFallback.Add(ghInt.Value);
            }

            // Flatten length fallback list
            List<int> lengthFallback = new List<int>();
            foreach (var p in lengthTree.Paths)
            {
                foreach (var ghInt in lengthTree[p])
                    lengthFallback.Add(ghInt.Value);
            }

            // Output trees
            var panelPlanesTree = new GH_Structure<GH_Plane>();        // Intrados (per vault -> per cell)
            var divisionPlanesTree = new GH_Structure<GH_Plane>();     // Each branch {branch; s; l} holds 4 planes per cell
            var transversalPlanesTree = new GH_Structure<GH_Plane>(); // Per vault branch -> transversal planes indexed by l
            var longitudinalPlanesTree = new GH_Structure<GH_Plane>(); // Per vault branch -> longitudinal planes indexed by s

            // For any global fallback lists (if needed), keep these
            // but we'll return trees to keep structure consistent
            foreach (GH_Path branchPath in vaultTree.Paths)
            {
                // Default values
                int spanDiv = 12;
                int lengthDiv = 8;

                // 1) Exact path match (preferred)
                if (spanTree.PathExists(branchPath) && spanTree[branchPath].Count > 0)
                    spanDiv = spanTree[branchPath][0].Value;

                // 2) Else try to match by branch index (last index of the path)
                else
                {
                    int branchIndex = (branchPath.Indices.Length > 0) ? branchPath.Indices[branchPath.Indices.Length - 1] : 0;
                    if (spanFallback.Count > branchIndex)
                        spanDiv = spanFallback[branchIndex];
                    else if (spanFallback.Count > 0)
                        spanDiv = spanFallback[0]; // broadcast first value
                }

                // Same for length
                if (lengthTree.PathExists(branchPath) && lengthTree[branchPath].Count > 0)
                    lengthDiv = lengthTree[branchPath][0].Value;
                else
                {
                    int branchIndex = (branchPath.Indices.Length > 0) ? branchPath.Indices[branchPath.Indices.Length - 1] : 0;
                    if (lengthFallback.Count > branchIndex)
                        lengthDiv = lengthFallback[branchIndex];
                    else if (lengthFallback.Count > 0)
                        lengthDiv = lengthFallback[0];
                }

                // Clamp values
                if (spanDiv < 3) spanDiv = 3;
                if (lengthDiv < 1) lengthDiv = 1;


                // collect surfaces in this branch
                List<Surface> surfaces = new List<Surface>();

                foreach (var ghSrf in vaultTree[branchPath])
                {
                    if (ghSrf == null || ghSrf.Value == null)
                        continue;

                    GeometryBase geo = ghSrf.Value;

                    if (geo is Surface srf)
                    {
                        // Already a surface — OK
                        surfaces.Add(srf);
                    }
                    else if (geo is Brep brep)
                    {
                        // If the brep has exactly 1 face, use its underlying surface
                        if (brep.Faces.Count == 1)
                        {
                            try
                            {
                                surfaces.Add(brep.Faces[0].UnderlyingSurface());
                            }
                            catch { }
                        }
                        else
                        {
                            // If multiple faces, try using the outermost face or skip
                            try
                            {
                                surfaces.Add(brep.Faces[0].UnderlyingSurface());
                            }
                            catch { }
                        }
                    }
                }

                if (surfaces.Count == 0) continue;

                // ---------- Build a single Brep/Surface for this branch ----------
                Brep loftBrep = null;
                try
                {
                    if (surfaces.Count == 1)
                    {
                        loftBrep = surfaces[0].ToBrep();
                    }
                    else
                    {
                        // Try gather profile curves from the surfaces' outer edges and loft them
                        var profileCurves = new List<Curve>();
                        foreach (var s in surfaces)
                        {
                            // Prefer the longest edge from each surface as profile
                            var edges = s.ToBrep().DuplicateEdgeCurves(true);
                            if (edges != null && edges.Length > 0)
                            {
                                // pick longest curve from this surface
                                var longest = edges.OrderByDescending(c => c.GetLength()).First();
                                profileCurves.Add(longest);
                            }
                        }

                        if (profileCurves.Count >= 2)
                        {
                            var loftBreps = Brep.CreateFromLoft(profileCurves, Point3d.Unset, Point3d.Unset, LoftType.Normal, false);
                            if (loftBreps != null && loftBreps.Length > 0)
                                loftBrep = loftBreps[0];
                        }
                    }
                }
                catch
                {
                    // if anything fails, skip this branch
                    continue;
                }

                if (loftBrep == null || loftBrep.Surfaces.Count == 0) continue;

                // Align normal to world Z (keeps original behaviour)
                Surface loftSrf = Components.PolylineUtils.AlignNormalToWorldZ(loftBrep.Surfaces[0]);

                // ---------- Extract boundary curves (arcs) ----------
                Curve[] boundaryCurvesArray = loftBrep.DuplicateEdgeCurves(true);
                List<Curve> boundaryCurves = (boundaryCurvesArray != null) ? boundaryCurvesArray.ToList() : new List<Curve>();
                if (boundaryCurves.Count < 2) continue;

                // Safer arc detection: take the two longest boundary curves (usually the two vault arcs)
                var arcCandidates = boundaryCurves.ToList();
                Curve arc1 = null;
                Curve arc2 = null;
                // find two sufficiently distinct longest curves
                if (arcCandidates.Count >= 2)
                {
                    arc1 = arcCandidates[1];
                    arc2 = arcCandidates[3];
                }
                if (arc1 == null || arc2 == null) continue;
                if (!arc1.IsValid || !arc2.IsValid) continue;

                // Orient arcs if you have a utility for that (keeps previous behaviour)
                try
                {
                    Components.Utils.OrientArcs(new List<Curve> { arc1, arc2 });
                }
                catch
                {
                    // ignore if utility fails
                }

                // Divide arcs spanwise
                List<Point3d> arc1_pts = new List<Point3d>();
                List<Point3d> arc2_pts = new List<Point3d>();
                double[] params1 = null;
                double[] params2 = null;

                try
                {
                    params1 = arc1.DivideByCount(spanDiv, true);
                    params2 = arc2.DivideByCount(spanDiv, true);
                    foreach (var t in params1) arc1_pts.Add(arc1.PointAt(t));
                    foreach (var t in params2) arc2_pts.Add(arc2.PointAt(t));
                }
                catch
                {
                    continue; // if division fails, skip branch
                }

                // Ensure we have the expected number of points
                if (arc1_pts.Count != arc2_pts.Count) continue;
                int spanPointCount = arc1_pts.Count; // expected spanDiv + 1

                // ---------- Build grid points [s][l] ----------
                List<List<Point3d>> grid = new List<List<Point3d>>();
                for (int s = 0; s < spanPointCount; s++)
                {
                    var row = new List<Point3d>();
                    // For each span position, line between corresponding points on arcs
                    Line spanLine = new Line(arc1_pts[s], arc2_pts[s]);
                    for (int l = 0; l <= lengthDiv; l++)
                    {
                        double factor = (lengthDiv == 0) ? 0.0 : (double)l / lengthDiv;
                        row.Add(spanLine.PointAt(factor));
                    }
                    grid.Add(row);
                }

                // ---------- Longitudinal planes (per s) ----------
                var longitudinalPlanes = new List<Plane>();
                double offsetDistance = 1.0; // as in original code
                for (int s = 0; s < spanPointCount; s++)
                {
                    // average across length for this span line
                    double sumX = 0, sumY = 0, sumZ = 0;
                    for (int l = 0; l <= lengthDiv; l++)
                    {
                        Point3d pt = grid[s][l];
                        sumX += pt.X; sumY += pt.Y; sumZ += pt.Z;
                    }
                    Point3d avgPt = new Point3d(sumX / (lengthDiv + 1), sumY / (lengthDiv + 1), sumZ / (lengthDiv + 1));

                    Point3d lastPt = grid[s][lengthDiv];

                    List<Point3d> movedPts = new List<Point3d>();
                    for (int l = 0; l <= lengthDiv; l++)
                    {
                        Point3d pt = grid[s][l];
                        double u, v;
                        if (loftSrf.ClosestPoint(pt, out u, out v))
                        {
                            Vector3d normal = loftSrf.NormalAt(u, v);
                            normal.Unitize();
                            Point3d ptMoved = pt + normal * offsetDistance;
                            movedPts.Add(ptMoved);
                        }
                        else
                        {
                            movedPts.Add(pt); // fallback
                        }
                    }

                    // average moved pts
                    double sumMX = 0, sumMY = 0, sumMZ = 0;
                    foreach (var mpt in movedPts)
                    {
                        sumMX += mpt.X; sumMY += mpt.Y; sumMZ += mpt.Z;
                    }
                    Point3d avgMovedPt = new Point3d(sumMX / movedPts.Count, sumMY / movedPts.Count, sumMZ / movedPts.Count);

                    Plane longitudinalPlane = new Plane(avgPt, lastPt, avgMovedPt);
                    longitudinalPlanes.Add(longitudinalPlane);

                    // Append to longitudinalPlanesTree under branch path
                    // Use path {branch; s} so they are grouped by branch and indexed by s
                    GH_Path outLp = new GH_Path(branchPath.Indices).AppendElement(s);
                    longitudinalPlanesTree.Append(new GH_Plane(longitudinalPlane), outLp);
                }

                // ---------- Transversal planes (per l) ----------
                var transversalPlanes = new List<Plane>();
                for (int l = 0; l <= lengthDiv; l++)
                {
                    List<Point3d> pts = new List<Point3d>();
                    List<Point3d> movedPts = new List<Point3d>();
                    for (int s = 0; s < spanPointCount; s++)
                    {
                        Point3d pt = grid[s][l];
                        pts.Add(pt);
                        double u, v;
                        if (loftSrf.ClosestPoint(pt, out u, out v))
                        {
                            // Original used Z axis offset for moved pts in transversal
                            Point3d ptMoved = pt + Vector3d.ZAxis * offsetDistance;
                            movedPts.Add(ptMoved);
                        }
                        else
                        {
                            movedPts.Add(pt);
                        }
                    }

                    // average original pts
                    double sumX = 0, sumY = 0, sumZ = 0;
                    foreach (var pt in pts)
                    {
                        sumX += pt.X; sumY += pt.Y; sumZ += pt.Z;
                    }
                    Point3d avgPt = new Point3d(sumX / pts.Count, sumY / pts.Count, sumZ / pts.Count);

                    Point3d lastPt = grid[spanPointCount - 1][l];

                    double sumMX = 0, sumMY = 0, sumMZ = 0;
                    foreach (var mpt in movedPts)
                    {
                        sumMX += mpt.X; sumMY += mpt.Y; sumMZ += mpt.Z;
                    }
                    Point3d avgMovedPt = new Point3d(sumMX / movedPts.Count, sumMY / movedPts.Count, sumMZ / movedPts.Count);

                    Plane transversalPlane = new Plane(avgPt, lastPt, avgMovedPt);
                    transversalPlanes.Add(transversalPlane);

                    // Append to transversalPlanesTree under branch path
                    // Use path {branch; l}
                    GH_Path outTp = new GH_Path(branchPath.Indices).AppendElement(l);
                    transversalPlanesTree.Append(new GH_Plane(transversalPlane), outTp);
                }

                // ---------- Panel (intrados) planes per cell ----------
                // For each cell s x l compute fit plane and then orient using loftSrf normal
                for (int s = 0; s < spanPointCount - 1; s++)
                {
                    for (int l = 0; l < lengthDiv; l++)
                    {
                        var a = grid[s][l];
                        var b = grid[s + 1][l];
                        var c = grid[s + 1][l + 1];
                        var d = grid[s][l + 1];
                        var pts = new Point3d[] { a, b, c, d };

                        Plane fitPlane;
                        Plane.FitPlaneToPoints(pts, out fitPlane);

                        Point3d avgPt = new Point3d(
                            (a.X + b.X + c.X + d.X) / 4.0,
                            (a.Y + b.Y + c.Y + d.Y) / 4.0,
                            (a.Z + b.Z + c.Z + d.Z) / 4.0
                        );

                        // surface normal fallback
                        Vector3d surfNormal = fitPlane.ZAxis;
                        double u, v;
                        if (loftSrf.ClosestPoint(avgPt, out u, out v))
                        {
                            try
                            {
                                surfNormal = loftSrf.NormalAt(u, v);
                                surfNormal.Unitize();
                            }
                            catch { }
                        }

                        Vector3d panelNormal = fitPlane.ZAxis;
                        panelNormal.Unitize();

                        // Flip the plane if its Z axis is opposed to the surface normal
                        if (!panelNormal.IsZero && !surfNormal.IsZero && panelNormal * surfNormal < 0)
                        {
                            fitPlane.Flip();
                            panelNormal = fitPlane.ZAxis;
                            panelNormal.Unitize();
                        }

                        // Try to align rotation similarly to the original algorithm
                        double angle = Vector3d.VectorAngle(fitPlane.XAxis, Plane.WorldXY.XAxis, fitPlane);
                        if (angle > Rhino.RhinoMath.ZeroTolerance)
                        {
                            fitPlane.Rotate(angle, panelNormal);
                        }

                        Plane panelPlane = new Plane(avgPt, fitPlane.XAxis, fitPlane.YAxis);

                        // Append to panelPlanesTree under branch path, use {branch; s; l}
                        GH_Path panelPath = new GH_Path(branchPath.Indices).AppendElement(s).AppendElement(l);
                        panelPlanesTree.Append(new GH_Plane(panelPlane), panelPath);

                        // Build the division planes list for this cell (start/end longitudinal & transversal)
                        // Use same branch path {branch; s; l} and append four planes
                        GH_Path divPath = new GH_Path(branchPath.Indices).AppendElement(s).AppendElement(l);
                        // safety: longitudinalPlanes has spanPointCount elements
                        divisionPlanesTree.Append(new GH_Plane(longitudinalPlanes[s]), divPath);       // Start longitudinal
                        divisionPlanesTree.Append(new GH_Plane(longitudinalPlanes[s + 1]), divPath);   // End longitudinal
                        divisionPlanesTree.Append(new GH_Plane(transversalPlanes[l]), divPath);        // Start transversal
                        divisionPlanesTree.Append(new GH_Plane(transversalPlanes[l + 1]), divPath);    // End transversal
                    }
                }
            } // end foreach branch

            panelPlanesTree.Simplify(GH_SimplificationMode.CollapseAllOverlaps);
            divisionPlanesTree.Simplify(GH_SimplificationMode.CollapseAllOverlaps);
            transversalPlanesTree.Simplify(GH_SimplificationMode.CollapseAllOverlaps);
            longitudinalPlanesTree.Simplify(GH_SimplificationMode.CollapseAllOverlaps);
            var tplanes = TreeUtils.TrimTreeDepth(transversalPlanesTree);
            var lplanes = TreeUtils.TrimTreeDepth(longitudinalPlanesTree);
            // Set outputs
            DA.SetDataTree(0, panelPlanesTree);
            DA.SetDataTree(1, divisionPlanesTree);
            DA.SetDataTree(2, tplanes);
            DA.SetDataTree(3, lplanes);
        }
    }
}
