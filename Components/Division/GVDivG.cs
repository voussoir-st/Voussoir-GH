using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
using Rhino.Geometry;
using Rhino.Geometry.Intersect;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

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

            pManager.AddNumberParameter(
                "V Divisions", "V2",
                "Number of voussoir divisions along the V direction (Profile direction)\nMinimum: 3",
                GH_ParamAccess.tree, 0.3);
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
            pManager.AddGenericParameter("V Planes", "Pv", "Division Planes with constant V-value.", GH_ParamAccess.list);
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
            GH_Structure<GH_Number> VTree2;

            DA.GetDataTree(1, out UTree);
            DA.GetDataTree(2, out VTree1);
            DA.GetDataTree(3, out VTree2);

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
                double VDiv2 = 0.3;

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

                //Debug.WriteLine($"surfaces: {breps.Count}");
                int n = 0;
                List<Curve> arcs = new List<Curve>();
                List<Curve> springerLines = new List<Curve>();
                List<Curve> notCandidates = new List<Curve>();

                foreach (var s in breps)
                {
                    n++;
                    //Debug.WriteLine($"\nsurface: {n}");

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


                    //Debug.WriteLine($"srfEdges: {srfEdges.ToList().Count}");

                    var lowest4 = vertices
                        .OrderBy(pt => pt.Z)
                        .Take(4)
                        .ToList();

                    Plane basePl;
                    Plane.FitPlaneToPoints(lowest4, out basePl);
                    //Debug.WriteLine($"dotprod: {basePl.ZAxis * Plane.WorldXY.ZAxis}");
                    //_previewBrep.Add(basePl);
                    if (basePl.ZAxis * Plane.WorldXY.ZAxis < 0)
                    {
                        //Debug.WriteLine($"dotprod: {basePl.ZAxis * Plane.WorldXY.ZAxis}");
                        basePl.Flip();
                    }

                    Curve arc = null;
                    List<Curve> spls = new List<Curve>();
                    int k = 0;
                    foreach (var edge in srfEdges)
                    {
                        k++;
                        //Debug.WriteLine($"edge: {k}");

                        // Check if both endpoints lie on the plane
                        bool endpointsOnPlane =
                            basePl.DistanceTo(edge.PointAtStart) <= tol &&
                            basePl.DistanceTo(edge.PointAtEnd) <= tol;
                        //Debug.WriteLine($"endpointsOnPlane: {endpointsOnPlane}");

                        if (!endpointsOnPlane)
                        {
                            // Store the non‑candidate edge
                            notCandidates.Add(edge.ToNurbsCurve());
                            continue;
                        }

                        // Now it's safe to treat this edge as a candidate
                        double midDist = basePl.DistanceTo(edge.PointAtNormalizedLength(0.5));
                        //Debug.WriteLine($"midDist: {midDist}");

                        if (midDist > tol)
                        {
                            arc = edge;                            
                        }
                        else
                        {
                            springerLines.Add(edge.ToNurbsCurve());
                        }
                    }

                    if (arc != null)
                    {
                        //_previewBrep.Add(arc.ToNurbsCurve());
                        //Debug.WriteLine($"arc: true");
                    }
                    //else Debug.WriteLine($"arc: false");

                    arcs.Add(arc);
                    //springerLines.AddRange(spls);
                }

                var arcLongest = arcs.OrderByDescending(l => l.GetLength()).First();

                Point3d[] UintPoints;
                arcLongest.DivideByCount(UDiv, true, out UintPoints);

                //intPoints.Select(p => new Point(p))
                //    .ToList()
                //    .ForEach(pt => _previewBrep.Add(pt));

                var intPointsHalf = UintPoints.ToList().Take(UintPoints.Length / 2);

                //intPointsHalf.Select(p => new Point(p))
                //   .ToList()
                //   .ForEach(pt => _previewBrep.Add(pt));

                List<Plane> splitPlanes = new List<Plane>();
                //_previewBrep.AddRange(springerLines);
                foreach (var p in intPointsHalf)
                {
                    var splitPlane = new Plane(p, Vector3d.ZAxis);
                    splitPlanes.Add(splitPlane);
                }
                //splitPlanes.RemoveAt(0);
                //_previewBrep.AddRange(splitPlanes);
                int i = 0;
                foreach (var s in breps)
                {
                    //Debug.WriteLine($"\nvault: {i}");
                    Curve arc = arcs[i];
                    List<Point3d> arcUpoints = new List<Point3d>();
                    List<Point3d> ncUpoints = new List<Point3d>();
                    var notCand = new List<Curve> { notCandidates[i * 2], notCandidates[i * 2 + 1] };
                    notCand = notCand
                        .OrderBy(u =>
                        {
                            double t;
                            arc.ClosestPoint(u.PointAtStart, out t);   // <-- THIS is the key difference
                            return t;
                        })
                        .ToList();

                    foreach (var p in splitPlanes)
                    {
                        var events = Intersection.CurvePlane(arc, p, tol);

                        List<Point3d> intersectionPoints = new List<Point3d>();

                        if (events != null)
                        {
                            foreach (var ev in events)
                                intersectionPoints.Add(ev.PointA);
                        }

                        arcUpoints.AddRange(intersectionPoints);

                        foreach (var v in notCand)
                        {
                            var evts = Intersection.CurvePlane(v, p, tol);

                            if (evts != null)
                            {
                                foreach (var ev in evts)
                                    ncUpoints.Add(ev.PointA);
                            }
                        }
                    }

                    ncUpoints = ncUpoints
                        .OrderBy(u =>
                        {
                            double t;
                            arc.ClosestPoint(u, out t);
                            return t;
                        })
                        .ToList();

                    arcUpoints = arcUpoints
                         .OrderBy(u =>
                         {
                             double t;
                             arc.ClosestPoint(u, out t);   // <-- THIS is the key difference
                             return t;
                         })
                         .ToList();
                    List<Curve> ULines = new List<Curve>();
                    for (int a = 0; a < ncUpoints.Count; a++)
                    {
                        var l = new Line(arcUpoints[a], ncUpoints[a]);
                        ULines.Add(new LineCurve(l));
                    }

                    var CrvLongest = ULines.OrderByDescending(l => l.GetLength()).First();

                    Point3d[] V1intPoints;

                    ULines = ULines
                         .OrderBy(u =>
                         {
                             double t;
                             arc.ClosestPoint(u.PointAtStart, out t);   // <-- THIS is the key difference
                             return t;
                         })
                         .ToList();

                    //interUCrvs.Insert(interUCrvs.Count, springerLines[i * 2]);
                    //interUCrvs.Insert(0, springerLines[i * 2 + 1]);

                    //_previewBrep.AddRange(ULines);

                    var p0 = arc.PointAtStart;
                    var p1 = arc.PointAtEnd;
                    var pmid = arc.PointAtNormalizedLength(0.5);

                    var avg = new Point3d(
                        (p0.X + p1.X + pmid.X) / 3.0,
                        (p0.Y + p1.Y + pmid.Y) / 3.0,
                        (p0.Z + p1.Z + pmid.Z) / 3.0
                    );

                    Plane Vplane = new Plane(avg, arc.PointAtStart, arc.PointAtEnd);

                    if (- tol > Vplane.DistanceTo(CrvLongest.PointAtStart) || Vplane.DistanceTo(CrvLongest.PointAtStart) > tol) CrvLongest.Reverse();

                    CrvLongest.DivideByCount(VDiv1, true, out V1intPoints);

                    //_previewBrep.Add(CrvLongest);

                    List<Plane> Vplanes = new List<Plane>();
                    var NC = Curve.JoinCurves(notCandidates);
                    foreach (var p in V1intPoints)
                    {
                        Plane vp = new Plane(p, Vplane.ZAxis);
                        Vplanes.Add(vp);
                        //_previewBrep.Add(vp);
                    }

                    for (int j = 0; j < ULines.Count - 1; j++)
                    {
                        //Debug.WriteLine($"\nline: {j}");

                        List<Point3d> intUPoints1 = new List<Point3d>();
                        List<Point3d> intUPoints2 = new List<Point3d>();
                        var c1 = ULines[j];
                        var c2 = ULines[j + 1];

                        foreach (var l in Vplanes)
                        {
                            var ev0 = Intersection.CurvePlane(c1, l, tol);
                            if (ev0 != null && ev0.Count > 0)
                                intUPoints1.Add(ev0[0].PointA);

                            var ev1 = Intersection.CurvePlane(c2, l, tol);
                            if (ev1 != null && ev1.Count > 0)
                                intUPoints2.Add(ev1[0].PointA);
                        }
                        
                        intUPoints1.Add(c1.PointAtEnd);
                        intUPoints2.Add(c2.PointAtEnd);
                        //Debug.WriteLine($"intUPoints1: {intUPoints1.Count}");
                        //Debug.WriteLine($"intUPoints2: {intUPoints2.Count}");
                        //Intrados
                        //_previewBrep.AddRange(intUPoints1.Select(p => new Point(p)));
                        //_previewBrep.AddRange(intUPoints2.Select(p => new Point(p)));

                        if (intUPoints1.Count < intUPoints2.Count)
                        {
                            int index1 = intUPoints1.Count - 2;
                            int index2 = intUPoints2.Count - 1;

                            // Remove all points between index1 and index2 (exclusive)
                            int start = index1 + 1;
                            int count = index2 - index1 - 1;

                            if (count > 0)
                                intUPoints2.RemoveRange(start, count);
                        }
                        else if (intUPoints1.Count > intUPoints2.Count)
                        {
                            int index1 = intUPoints1.Count - 1;
                            int index2 = intUPoints2.Count - 2;

                            // Remove all points between index1 and index2 (exclusive)
                            int start = index2 + 1;
                            int count = index1 - index2 - 1;

                            if (count > 0)
                                intUPoints1.RemoveRange(start, count);
                        }

                        //_previewBrep.AddRange(intUPoints1.Select(p => new Point(p)));
                        //_previewBrep.AddRange(intUPoints2.Select(p => new Point(p)));
                        List<Plane> intradosPlanes = new List<Plane>();
                        
                        for (int k = 0; k < intUPoints1.Count - 1; k++)
                        {
                            var ip1 = intUPoints1[k];
                            var ip2 = intUPoints1[k + 1];
                            var ip3 = intUPoints2[k + 1];
                            var ip4 = intUPoints2[k];

                            List<Point3d> ipPoints1 = new List<Point3d> { ip1, ip2, ip3, ip4 , ip1};

                            var ipavg = new Point3d(
                                (ip1.X + ip2.X + ip3.X + ip4.X) / 4.0,
                                (ip1.Y + ip2.Y + ip3.Y + ip4.Y) / 4.0,
                                (ip1.Z + ip2.Z + ip3.Z + ip4.Z) / 4.0
                            );

                            Plane iPlane;
                            Plane.FitPlaneToPoints(ipPoints1, out iPlane);
                            iPlane = new Plane(ipavg, iPlane.ZAxis);
                            intradosPlanes.Add(iPlane);

                            // Midpoints
                            var ip1_2 = new Point3d((ip1.X + ip2.X) / 2, (ip1.Y + ip2.Y) / 2, (ip1.Z + ip2.Z) / 2);
                            var ip2_3 = new Point3d((ip2.X + ip3.X) / 2, (ip2.Y + ip3.Y) / 2, (ip2.Z + ip3.Z) / 2);
                            var ip3_4 = new Point3d((ip3.X + ip4.X) / 2, (ip3.Y + ip4.Y) / 2, (ip3.Z + ip4.Z) / 2);
                            var ip4_1 = new Point3d((ip4.X + ip1.X) / 2, (ip4.Y + ip1.Y) / 2, (ip4.Z + ip1.Z) / 2);

                            if (s.ClosestPoint(ip1_2, out Point3d aa, out ComponentIndex ab, out double ac, out double ad, 0.0, out Vector3d normal_ip1_2))
                            {
                                normal_ip1_2.Unitize();

                                if (normal_ip1_2 * iPlane.ZAxis < 0)
                                    normal_ip1_2.Reverse();
                            }
                            if (s.ClosestPoint(ip2_3, out Point3d ba, out ComponentIndex bb, out double bc, out double bd, 0.0, out Vector3d normal_ip2_3))
                            {
                                normal_ip2_3.Unitize();

                                if (normal_ip2_3 * iPlane.ZAxis < 0)
                                    normal_ip2_3.Reverse();
                            }
                            if (s.ClosestPoint(ip3_4, out Point3d ca, out ComponentIndex cb, out double cc, out double cd, 0.0, out Vector3d normal_ip3_4))
                            {
                                normal_ip3_4.Unitize();

                                if (normal_ip3_4 * iPlane.ZAxis < 0)
                                    normal_ip3_4.Reverse();
                            }
                            if (s.ClosestPoint(ip4_1, out Point3d da, out ComponentIndex db, out double dc, out double dd, 0.0, out Vector3d normal_ip4_1))
                            {
                                normal_ip4_1.Unitize();

                                if (normal_ip4_1 * iPlane.ZAxis < 0)
                                    normal_ip4_1.Reverse();
                            }

                            Vector3d n1 = (normal_ip4_1 + normal_ip1_2) * 0.5;
                            n1.Unitize();
                            Vector3d n2 = (normal_ip1_2 + normal_ip2_3) * 0.5;
                            Vector3d n3 = (normal_ip2_3 + normal_ip3_4) * 0.5;
                            Vector3d n4 = (normal_ip3_4 + normal_ip4_1) * 0.5;

                            n2.Unitize();
                            n3.Unitize();
                            n4.Unitize();

                            var ip5 = ip1 + n1 * VDiv2;
                            var ip6 = ip2 + n2 * VDiv2;
                            var ip7 = ip3 + n3 * VDiv2;
                            var ip8 = ip4 + n4 * VDiv2;

                            List<Point3d> ipPoints2 = new List<Point3d> { ip5, ip6, ip7, ip8, ip5 };

                            var intPl1 = new PolylineCurve(ipPoints1);
                            var intPl2 = new PolylineCurve(ipPoints2);
                            _previewBrep.Add(intPl1);
                            //_previewBrep.Add(intPl2);

                        }
                    }
                    i++;
                }
            }
            DA.SetDataList(5, _previewBrep);
        }
    }
}
