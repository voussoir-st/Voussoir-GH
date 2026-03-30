using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
using Rhino;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Net;
using VoussoirPlugin03.Components.Springers;
using VoussoirPlugin03.Properties;

namespace VoussoirPlugin03.Components.Division
{
    public class GroinVaultDivisionComponent : GH_Component
    {
        public GroinVaultDivisionComponent()
          : base(
                "Groin Vault Division - Grid",
                "GVDivG",
                "Divides a vault defined by two arcs into spanwise and lengthwise voussoirs. Accepts a tree of surfaces (each branch one or more surfaces).",
                "Voussoir",
                "Division"
                )
        { }

        public override Guid ComponentGuid => new Guid("BC98E9F2-CD3B-4C41-ADFF-FD189794437C");

        //protected override System.Drawing.Bitmap Icon
        //{
        //    get
        //    {
        //        return VoussoirPlugin03.Properties.Resources.VDiv;
        //    }
        //}

        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            // Now accepts a tree of surfaces
            pManager.AddSurfaceParameter("Vault Surface", "S", "Surface(s) that define the vault. Provide as a tree where each branch corresponds to one vault (1 or multiple surfaces).", GH_ParamAccess.tree);
            pManager.AddIntegerParameter(
                "U Divisions", "U",
                "Number of voussoir divisions along the U direction (Springer Lines direction)\nMinimum: 1",
                GH_ParamAccess.tree, 8);

            pManager.AddIntegerParameter(
                "V Divisions", "V1",
                "Number of voussoir divisions along the V direction (Profile direction)\nMinimum: 3",
                GH_ParamAccess.tree, 12);

            pManager.AddIntegerParameter(
                "V Divisions", "V2",
                "Number of voussoir divisions along the V direction (Profile direction)\nMinimum: 3",
                GH_ParamAccess.tree, 12);
        }

        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddPlaneParameter("Intrados Planes", "Pi", "Intrados planar vault panels (tree matching input branches).", GH_ParamAccess.tree);
            pManager.AddPlaneParameter("Division planes", "Pd", "Planes of each Voussoir Contact Surface (tree matching input branches).", GH_ParamAccess.tree);
            pManager.HideParameter(1);
            pManager.AddTextParameter("Boundaries", "B", "Voussoir Boundaries (Indexes of Division Planes).", GH_ParamAccess.tree);
            pManager.HideParameter(2);
            pManager.AddPlaneParameter("U Planes", "Pu", "Division Planes with constant U-value.", GH_ParamAccess.list);
            pManager.HideParameter(3);
            pManager.AddPlaneParameter("V Planes", "Pv", "Division Planes with constant V-value.", GH_ParamAccess.list);
            pManager.HideParameter(4);
            pManager.AddGeometryParameter("V Planes", "Pv", "Division Planes with constant V-value.", GH_ParamAccess.list);
        }
        private List<Point3d> _previewPts = new List<Point3d>();
        private List<Curve> _previewCrvs = new List<Curve>();
        //private List<GeometryBase> _previewBrep = new List<GeometryBase>();
        public override bool IsPreviewCapable => true;
        public override void DrawViewportWires(IGH_PreviewArgs args)
        {
            base.DrawViewportWires(args);

            // Draw points
            foreach (var pt in _previewPts)
            {
                args.Display.DrawPoint(pt, Rhino.Display.PointStyle.RoundSimple, 3, System.Drawing.Color.DarkRed);
            }

            // Draw curves
            foreach (var crv in 
                _previewCrvs)
            {
                if (crv != null && crv.IsValid)
                    args.Display.DrawCurve(crv, System.Drawing.Color.DarkRed, 2);
            }
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            _previewPts.Clear();
            _previewCrvs.Clear();
            List<GeometryBase> _previewBrep = new List<GeometryBase>();

            double tol = 0.0001;
            var vaultTree = new GH_Structure<GH_Surface>();
            if (!DA.GetDataTree(0, out vaultTree)) return;
            // Read span and length divisions as trees
            GH_Structure<GH_Integer> UTree;
            GH_Structure<GH_Integer> VTree1;
            GH_Structure<GH_Integer> VTree2;

            DA.GetDataTree(2, out UTree);
            DA.GetDataTree(1, out VTree1);
            DA.GetDataTree(1, out VTree2);

            // Flatten span fallback list
            List<int> UFallback = new List<int>();
            foreach (var p in UTree.Paths)
            {
                foreach (var ghInt in UTree[p])
                    UFallback.Add(ghInt.Value);
            }

            // Flatten length fallback list
            List<int> V1Fallback = new List<int>();
            foreach (var p in VTree1.Paths)
            {
                foreach (var ghInt in VTree1[p])
                    V1Fallback.Add(ghInt.Value);
            }
            List<int> V2Fallback = new List<int>();
            foreach (var p in VTree1.Paths)
            {
                foreach (var ghInt in VTree1[p])
                    V2Fallback.Add(ghInt.Value);
            }

            foreach (GH_Path branchPath in vaultTree.Paths)
            {
                // Default values
                int UDiv = 12;
                int VDiv1 = 8;
                int VDiv2 = 8;

                // 1) Exact path match (preferred)
                if (UTree.PathExists(branchPath) && UTree[branchPath].Count > 0)
                    UDiv = UTree[branchPath][0].Value;

                // 2) Else try to match by branch index (last index of the path)
                else
                {
                    int branchIndex = (branchPath.Indices.Length > 0) ? branchPath.Indices[branchPath.Indices.Length - 1] : 0;
                    if (UFallback.Count > branchIndex)
                        UDiv = UFallback[branchIndex];
                    else if (UFallback.Count > 0)
                        UDiv = UFallback[0]; // broadcast first value
                }

                // Same for length
                if (VTree1.PathExists(branchPath) && VTree1[branchPath].Count > 0)
                    VDiv1 = VTree1[branchPath][0].Value;
                else
                {
                    int branchIndex = (branchPath.Indices.Length > 0) ? branchPath.Indices[branchPath.Indices.Length - 1] : 0;
                    if (V1Fallback.Count > branchIndex)
                        VDiv1 = V1Fallback[branchIndex];
                    else if (V1Fallback.Count > 0)
                        VDiv1 = V1Fallback[0];
                }
                // Same for length
                if (VTree2.PathExists(branchPath) && VTree2[branchPath].Count > 0)
                    VDiv2 = VTree2[branchPath][0].Value;
                else
                {
                    int branchIndex = (branchPath.Indices.Length > 0) ? branchPath.Indices[branchPath.Indices.Length - 1] : 0;
                    if (V2Fallback.Count > branchIndex)
                        VDiv2 = V2Fallback[branchIndex];
                    else if (V2Fallback.Count > 0)
                        VDiv2 = V2Fallback[0];
                }

                // Clamp values
                if (UDiv < 3) UDiv = 3;
                if (VDiv1 < 1) VDiv1 = 1;
                if (VDiv2 < 1) VDiv2 = 1;

                // collect surfaces in this branch
                List<Brep> breps = new List<Brep>();

                foreach (var ghSrf in vaultTree[branchPath])
                {
                    if (ghSrf == null || ghSrf.Value == null)
                        continue;

                    GeometryBase geo = ghSrf.Value;

                    Brep b = null;

                    if (geo is Brep brep)
                        b = brep;
                    else if (geo is Surface srf)
                        b = srf.ToBrep();
                    else if (geo is Extrusion ext)
                        b = ext.ToBrep();
                    else
                        continue;

                    if (b != null)
                        breps.Add(b);
                }

                Debug.WriteLine($"surfaces: {breps.Count}");

                foreach (var s in breps)
                {
                    
                    var srfEdges = s.Edges;
                    //_previewCrvs.AddRange(srfEdges);
                    //_previewBrep.AddRange(srfEdges);

                    // 2. Join them into a single curve (or as few as possible)
                    Curve[] joined = Curve.JoinCurves(srfEdges, tol);

                    if (joined == null || joined.Length == 0)
                        continue;

                    Curve boundary = joined[0];   // vault surfaces should give 1 closed curve

                    // 3. Extract vertices (start/end points)
                    // For closed curves, Start == End, so we get all segment endpoints
                    List<Point3d> vertices = new List<Point3d>();

                    if (boundary is PolyCurve poly)
                    {
                        foreach (Curve seg in poly.Explode())
                        {
                            vertices.Add(seg.PointAtStart);
                            vertices.Add(seg.PointAtEnd);
                        }
                    }
                    else
                    {
                        // Single segment curve
                        vertices.Add(boundary.PointAtStart);
                        vertices.Add(boundary.PointAtEnd);
                    }

                    // Optional: remove duplicates (closed curves repeat points)
                    vertices = vertices.Distinct().ToList();


                    Debug.WriteLine($"srfEdges: {srfEdges.ToList().Count}");

                    var lowest4 = vertices
                        .OrderBy(pt => pt.Z)
                        .Take(4)
                        .ToList();

                    Plane basePl;
                    Plane.FitPlaneToPoints(lowest4, out basePl);

                    Curve arc = null;

                    foreach (var edge in srfEdges)
                    {
                        // Check if both endpoints lie on the plane
                        bool endpointsOnPlane =
                            basePl.DistanceTo(edge.PointAtStart) <= tol &&
                            basePl.DistanceTo(edge.PointAtEnd) <= tol;

                        if (!endpointsOnPlane)
                            continue;

                        // Now it's safe to treat this edge as a candidate
                        double midDist = basePl.DistanceTo(edge.PointAt(0.5));

                        if (midDist > tol)
                        {
                            arc = edge;
                            break; // assuming only one arc edge exists
                        }
                    }

                    if (arc != null) _previewBrep.Add(arc.ToNurbsCurve());

                }


            }
            DA.SetDataList(5, _previewBrep);
        }
    }
}
