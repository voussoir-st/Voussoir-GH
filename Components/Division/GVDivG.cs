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

        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                return VoussoirPlugin03.Properties.Resources.GVDivG;
            }
        }

        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {            
            pManager.AddSurfaceParameter("Vault Surface", "S", "Surface(s) that define the vault. Provide as a tree where each branch corresponds to one vault (1 or multiple surfaces).", GH_ParamAccess.tree);
            
            pManager.AddIntegerParameter("U Divisions", "U1", "Number of voussoir divisions along the U direction (Springer Lines direction)\nMinimum: 1", GH_ParamAccess.tree, 8);
            
            pManager.AddIntegerParameter("U Divisions", "U2", "Number of voussoir divisions along the V direction (Profile direction)\nMinimum: 3", GH_ParamAccess.tree, 8);
            
            pManager.AddIntegerParameter("V Divisions", "V", "Number of voussoir divisions along the V direction (Profile direction)\nMinimum: 3", GH_ParamAccess.tree, 12);
            
            pManager.AddNumberParameter("Division Tolerance", "T", "Minimum distance from the Groin line to the second closest voussoir.", GH_ParamAccess.tree, 0.1);
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
            pManager.AddGenericParameter("V Planes", "Pv", "Division Planes with constant V-value.", GH_ParamAccess.list);
            pManager.HideParameter(4);
            //pManager.AddGenericParameter("V Planes", "Pv", "Division Planes with constant V-value.", GH_ParamAccess.list);
        }
        //private List<Point3d> _previewPts = new List<Point3d>();
        //private List<Curve> _previewCrvs = new List<Curve>();
        //private List<GeometryBase> _previewBrep = new List<GeometryBase>();
        public override bool IsPreviewCapable => true;
        //public override void DrawViewportWires(IGH_PreviewArgs args)
        //{
        //    base.DrawViewportWires(args);

        //    // Draw points
        //    foreach (var pt in _previewPts)
        //    {
        //        args.Display.DrawPoint(pt, Rhino.Display.PointStyle.RoundSimple, 3, System.Drawing.Color.DarkRed);
        //    }

        //    // Draw curves
        //    foreach (var crv in
        //        _previewCrvs)
        //    {
        //        if (crv != null && crv.IsValid)
        //            args.Display.DrawCurve(crv, System.Drawing.Color.DarkRed, 2);
        //    }
        //}

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // Output trees
            var panelPlanesTree = new GH_Structure<GH_Plane>();        // Intrados (per vault -> per cell)
            var divisionPlanesTree = new GH_Structure<GH_Plane>();     // Each branch {branch; s; l} holds 4 planes per cell
            var boundaries = new GH_Structure<GH_String>();
            var transversalPlanesTree = new GH_Structure<GH_Plane>(); // Per vault branch -> transversal planes indexed by l
            var longitudinalPlanesTree = new GH_Structure<GH_Plane>(); // Per vault branch -> longitudinal planes indexed by s

            //_previewPts.Clear();
            //_previewCrvs.Clear();
            List<Plane> _previewBrep = new List<Plane>();
            List<GeometryBase> _previewCrvs = new List<GeometryBase>();

            double tol = 0.0001;
            var vaultTree = new GH_Structure<GH_Surface>();
            if (!DA.GetDataTree(0, out vaultTree)) return;
            // Read span and length divisions as trees
            GH_Structure<GH_Integer> UTree1; 
            GH_Structure<GH_Integer> UTree2; 
            GH_Structure<GH_Integer> VTree;

            DA.GetDataTree(1, out UTree1);
            DA.GetDataTree(2, out UTree2);
            DA.GetDataTree(3, out VTree);
            GH_Structure<GH_Number> tolerances;
            DA.GetDataTree(4, out tolerances);

            List<double> tolFallback = new List<double>();
            foreach (var p in tolerances.Paths)
            {
                foreach (var ghInt in tolerances[p])
                    tolFallback.Add(ghInt.Value);
            }

            // Flatten span fallback list
            List<int> U1Fallback = new List<int>();
            foreach (var p in UTree1.Paths)
            {
                foreach (var ghInt in UTree1[p])
                    U1Fallback.Add(ghInt.Value);
            }

            // Flatten length fallback list
            List<int> U2Fallback = new List<int>();
            foreach (var p in UTree2.Paths)
            {
                foreach (var ghInt in UTree2[p])
                    U2Fallback.Add(ghInt.Value);
            }
            List<int> VFallback = new List<int>();
            foreach (var p in VTree.Paths)
            {
                foreach (var ghInt in VTree[p])
                    VFallback.Add(ghInt.Value);
            }

            int q = 0;
            foreach (GH_Path branchPath in vaultTree.Paths)
            {
                // Default values
                int U1Div = 12;
                int U2Div = 12;
                int VDiv = 8;
                double groinTol = 0.1;

                // 1) Exact path match (preferred)
                if (tolerances.PathExists(branchPath) && tolerances[branchPath].Count > 0)
                    groinTol = tolerances[branchPath][0].Value;

                // 2) Else try to match by branch index (last index of the path)
                else
                {
                    int branchIndex = (branchPath.Indices.Length > 0) ? branchPath.Indices[branchPath.Indices.Length - 1] : 0;
                    if (tolFallback.Count > branchIndex)
                        groinTol = tolFallback[branchIndex];
                    else if (tolFallback.Count > 0)
                        groinTol = tolFallback[0]; // broadcast first value
                }

                // 1) Exact path match (preferred)
                if (UTree1.PathExists(branchPath) && UTree1[branchPath].Count > 0)
                    U1Div = UTree1[branchPath][0].Value;

                // 2) Else try to match by branch index (last index of the path)
                else
                {
                    int branchIndex = (branchPath.Indices.Length > 0) ? branchPath.Indices[branchPath.Indices.Length - 1] : 0;
                    if (U1Fallback.Count > branchIndex)
                        U1Div = U1Fallback[branchIndex];
                    else if (U1Fallback.Count > 0)
                        U1Div = U1Fallback[0]; // broadcast first value
                }

                // Same for length
                if (UTree2.PathExists(branchPath) && UTree2[branchPath].Count > 0)
                    U2Div = UTree1[branchPath][0].Value;
                else
                {
                    int branchIndex = (branchPath.Indices.Length > 0) ? branchPath.Indices[branchPath.Indices.Length - 1] : 0;
                    if (U2Fallback.Count > branchIndex)
                        U2Div = U2Fallback[branchIndex];
                    else if (U2Fallback.Count > 0)
                        U2Div = U2Fallback[0];
                }
                // Same for length
                if (VTree.PathExists(branchPath) && VTree[branchPath].Count > 0)
                    VDiv = VTree[branchPath][0].Value;
                else
                {
                    int branchIndex = (branchPath.Indices.Length > 0) ? branchPath.Indices[branchPath.Indices.Length - 1] : 0;
                    if (VFallback.Count > branchIndex)
                        VDiv = VFallback[branchIndex];
                    else if (VFallback.Count > 0)
                        VDiv = VFallback[0];
                }

                // Clamp values
                if (U1Div < 1) U1Div = 1;
                if (U2Div < 1) U2Div = 1;
                if (VDiv < 3) VDiv = 3;

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
                arcLongest.DivideByCount(VDiv, true, out UintPoints);

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
                    var surfaceBoundaries = new List<String>();
                    var vaultdivisionPlanesTree = new List<Plane>();
                    var uPlanestree = new List<Plane>();
                    var vPlanestree = new List<Plane>();
                    var surf = s.Faces[0].UnderlyingSurface();
                    var sAlign = BaseSurface.PolylineUtils.AlignNormalToWorldZ(surf);
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
                    if (new Vector3d(notCand[0].PointAtEnd - notCand[0].PointAtStart) * Vector3d.ZAxis < 0) notCand[0].Reverse();

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
                        //_previewCrvs.AddRange(intersectionPoints.Select(a => new Point(a)));

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

                    var nc = Curve.JoinCurves(notCand)[0];
                    if (new Vector3d(nc.PointAtEnd - nc.PointAtStart) * new Vector3d(arc.PointAtEnd - arc.PointAtStart) < 0) 
                        nc.Reverse();

                    ncUpoints = ncUpoints
                        .OrderBy(u =>
                        {
                            double t;
                            nc.ClosestPoint(u, out t);
                            return t;
                        })
                        .ToList();

                    arcUpoints = arcUpoints
                         .OrderBy(u =>
                         {
                             double t;
                             arc.ClosestPoint(u, out t);
                             return t;
                         })
                         .ToList();

                    if (VDiv % 2 == 0)
                    {
                        arcUpoints.Insert(VDiv / 2, arc.PointAtNormalizedLength(0.5));
                        ncUpoints.Insert(VDiv / 2, notCand[0].PointAtEnd);
                    }

                    //_previewCrvs.AddRange(ncUpoints.Select(a => new Point(a)));

                    List<Curve> ULines = new List<Curve>();
                    for (int a = 0; a < ncUpoints.Count; a++)
                    {                        
                        ULines.Add(new LineCurve(new Line(arcUpoints[a], ncUpoints[a])));
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

                    //_previewCrvs.AddRange(ULines);

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

                    if (Params.Input[2].SourceCount > 0 && i == 2 || i== 3) CrvLongest.DivideByCount(U2Div, true, out V1intPoints);
                    else CrvLongest.DivideByCount(U1Div, true, out V1intPoints);

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
                        
                        intUPoints1 = BaseSurface.PolylineUtils.DistinctByTolerance(intUPoints1, groinTol);
                        intUPoints2 = BaseSurface.PolylineUtils.DistinctByTolerance(intUPoints2, groinTol);

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

                        //_previewCrvs.AddRange(intUPoints1.Select(p => new Point(p)));
                        //_previewCrvs.AddRange(intUPoints2.Select(p => new Point(p)));
                        List<Plane> intradosPlanes = new List<Plane>();
                        
                        for (int k = 0; k < intUPoints1.Count - 1; k++)
                        {
                            var ip1 = intUPoints1[k];
                            var ip2 = intUPoints1[k + 1];
                            var ip3 = intUPoints2[k + 1];
                            var ip4 = intUPoints2[k];

                            List<Point3d> ipPoints1 = new List<Point3d> { ip1, ip2, ip3, ip4, ip1 };

                            var ipavg = new Point3d(
                                (ip1.X + ip2.X + ip3.X + ip4.X) / 4.0,
                                (ip1.Y + ip2.Y + ip3.Y + ip4.Y) / 4.0,
                                (ip1.Z + ip2.Z + ip3.Z + ip4.Z) / 4.0
                            );

                            Plane iPlane;
                            Plane.FitPlaneToPoints(ipPoints1, out iPlane);

                            iPlane = new Plane(ipavg, iPlane.ZAxis);

                            double u, v;
                            sAlign.ClosestPoint(iPlane.Origin, out u, out v);
                            iPlane = new Plane(iPlane.Origin, sAlign.NormalAt(u, v));                           
                            
                            // Midpoints
                            var ip1_2 = new Point3d((ip1.X + ip2.X) / 2, (ip1.Y + ip2.Y) / 2, (ip1.Z + ip2.Z) / 2);
                            var ip2_3 = new Point3d((ip2.X + ip3.X) / 2, (ip2.Y + ip3.Y) / 2, (ip2.Z + ip3.Z) / 2);
                            var ip3_4 = new Point3d((ip3.X + ip4.X) / 2, (ip3.Y + ip4.Y) / 2, (ip3.Z + ip4.Z) / 2);
                            var ip4_1 = new Point3d((ip4.X + ip1.X) / 2, (ip4.Y + ip1.Y) / 2, (ip4.Z + ip1.Z) / 2);

                            var ip1N = Components.BaseSurface.PolylineUtils.GetBrepNormalAtPoint(s.Faces[0].UnderlyingSurface().ToBrep(), ip1, 1);
                            var ip2N = Components.BaseSurface.PolylineUtils.GetBrepNormalAtPoint(s.Faces[0].UnderlyingSurface().ToBrep(), ip2, 1);
                            var ip3N = Components.BaseSurface.PolylineUtils.GetBrepNormalAtPoint(s.Faces[0].UnderlyingSurface().ToBrep(), ip3, 1);
                            var ip4N = Components.BaseSurface.PolylineUtils.GetBrepNormalAtPoint(s.Faces[0].UnderlyingSurface().ToBrep(), ip4, 1);

                            Vector3d ip1_2N_avg = new Vector3d(
                                (ip1N.X + ip2N.X) / 2.0,
                                (ip1N.Y + ip2N.Y) / 2.0,
                                (ip1N.Z + ip2N.Z) / 2.0);
                            Vector3d ip2_3N_avg = new Vector3d(
                                (ip2N.X + ip3N.X) / 2.0,
                                (ip2N.Y + ip3N.Y) / 2.0,
                                (ip2N.Z + ip3N.Z) / 2.0);
                            Vector3d ip3_4N_avg = new Vector3d(
                                (ip3N.X + ip4N.X) / 2.0,
                                (ip3N.Y + ip4N.Y) / 2.0,
                                (ip3N.Z + ip4N.Z) / 2.0);
                            Vector3d ip4_1N_avg = new Vector3d(
                                (ip4N.X + ip1N.X) / 2.0,
                                (ip4N.Y + ip1N.Y) / 2.0,
                                (ip4N.Z + ip1N.Z) / 2.0);
                            Vector3d ipN_avg = new Vector3d(
                                (ip1N.X + ip2N.X + ip3N.X + ip4N.X) / 4.0,
                                (ip1N.Y + ip2N.Y + ip3N.Y + ip4N.Y) / 4.0,
                                (ip1N.Z + ip2N.Z + ip3N.Z + ip4N.Z) / 4.0);

                            iPlane = new Plane(iPlane.Origin, ipN_avg);
                            intradosPlanes.Add(iPlane);

                            var P1_2 = new Plane(ip1_2, new Vector3d(ip2 - ip1), ip1_2N_avg);
                            var P2_3 = new Plane(ip2_3, new Vector3d(ip3 - ip2), ip2_3N_avg);
                            var P3_4 = new Plane(ip3_4, new Vector3d(ip4 - ip3), ip3_4N_avg);
                            var P4_1 = new Plane(ip4_1, new Vector3d(ip1 - ip4), ip4_1N_avg);

                            var intPl1 = new PolylineCurve(ipPoints1);
                            if (k == intUPoints1.Count - 2 && j != ULines.Count / 2 - 1)
                            {
                                //Debug.WriteLine("true");

                                if (breps.Count < 2)
                                    return;

                                // --- Find the two closest breps to ip2 ---
                                var twoClosestToIp2 =
                                    breps
                                        .Select(b =>
                                        {
                                            Point3d cp = b.ClosestPoint(ip2);
                                            double d = ip2.DistanceTo(cp);
                                            return new { Brep = b, Dist = d };
                                        })
                                        .OrderBy(x => x.Dist)
                                        .Take(2)
                                        .Select(x => x.Brep)
                                        .ToList();

                                var b0 = twoClosestToIp2[0];
                                var b1 = twoClosestToIp2[1];

                                // --- Normals at ip2 ---
                                double na1, nb1;
                                b0.Faces[0].ClosestPoint(ip2, out na1, out nb1);
                                var n1 = b0.Faces[0].NormalAt(na1, nb1);
                                double na2, nb2;
                                b1.Faces[0].ClosestPoint(ip2, out na2, out nb2);
                                var n2 = b1.Faces[0].NormalAt(na2, nb2);

                                var l1 = new Line(ip2, ip2 + n1 * 0.3);
                                var l2 = new Line(ip2, ip2 + n2 * 0.3);

                                // --- Normals at ip3 ---
                                double na3, nb3;
                                b0.Faces[0].ClosestPoint(ip3, out na3, out nb3);
                                var n3 = b0.Faces[0].NormalAt(na3, nb3);
                                double na4, nb4;
                                b1.Faces[0].ClosestPoint(ip3, out na4, out nb4);
                                var n4 = b1.Faces[0].NormalAt(na4, nb4);

                                var l3 = new Line(ip3, ip3 + n3 * 0.3);
                                var l4 = new Line(ip3, ip3 + n4 * 0.3);

                                // --- Average all four normals ---
                                Vector3d nAvg = new Vector3d(
                                    (n1.X + n2.X + n3.X + n4.X) / 4.0,
                                    (n1.Y + n2.Y + n3.Y + n4.Y) / 4.0,
                                    (n1.Z + n2.Z + n3.Z + n4.Z) / 4.0
                                );

                                // --- Preview ---
                                //_previewCrvs.Add(new LineCurve(l1));
                                //_previewCrvs.Add(new LineCurve(l2));
                                //_previewCrvs.Add(new LineCurve(l3));
                                //_previewCrvs.Add(new LineCurve(l4));

                                // --- Final plane ---
                                P2_3 = new Plane(ip2_3, new Vector3d(ip3 - ip2), nAvg);
                            }

                            //var srfpl1 = intPl1.PullToBrepFace(s.Faces[0], tol);
                            _previewCrvs.Add(intPl1.ToNurbsCurve());
                            //_previewBrep.Add(iPlane);

                            panelPlanesTree.Append(new GH_Plane(iPlane), new GH_Path(q, i));

                            int count = vaultdivisionPlanesTree != null ? vaultdivisionPlanesTree.Count : 0;
                            boundaries.Append(new GH_String($"B{{{count};{count + 1};{count + 2};{count + 3}}}"), new GH_Path(q, i));

                            vaultdivisionPlanesTree.Add(P1_2);
                            uPlanestree.Add(P1_2);
                            vaultdivisionPlanesTree.Add(P2_3);
                            
                            vaultdivisionPlanesTree.Add(P3_4);
                            if (j == ULines.Count - 2) uPlanestree.Add(P3_4);
                            vaultdivisionPlanesTree.Add(P4_1);
                            vPlanestree.Add(P4_1);
                            if (k == intUPoints1.Count - 2) vPlanestree.Add(P2_3);
                        }
                    }
                    divisionPlanesTree.AppendRange(vaultdivisionPlanesTree.Select(p => new GH_Plane(p)), new GH_Path(q, i));
                    transversalPlanesTree.AppendRange(vPlanestree.Select(p => new GH_Plane(p)), new GH_Path(q, i));
                    longitudinalPlanesTree.AppendRange(uPlanestree.Select(p => new GH_Plane(p)), new GH_Path(q, i));
                    i++;


                }
                q++;
            }
            DA.SetDataTree(0, panelPlanesTree);
            DA.SetDataTree(1, divisionPlanesTree);
            DA.SetDataTree(2, boundaries);
            DA.SetDataTree(3, transversalPlanesTree);
            DA.SetDataList(4, longitudinalPlanesTree);
            //DA.SetDataList(5, _previewBrep);            
        }
    }
}
