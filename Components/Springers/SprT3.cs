using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
using Rhino;
using Rhino.Geometry;
using Rhino.Geometry.Intersect;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using VoussoirPlugin03.Components.BaseSurface;

namespace VoussoirPlugin03.Components.Springers
{
    public class VoussoirCreate1 : GH_Component
    {
        public VoussoirCreate1()
            : base("Springer - Trapezoid", "SprT",
                  "Create a simple trapezoid Springer",
                  "Voussoir", "Springer")
        { }

        public override Guid ComponentGuid => new Guid("FC88F9F2-CD3B-4C41-ADFF-FD189794137C");

        protected override System.Drawing.Bitmap Icon => VoussoirPlugin03.Properties.Resources.SpringerT;

        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddBrepParameter("Vault Surface", "S", "Surface that defines the vault", GH_ParamAccess.tree);
            pManager.AddCurveParameter("Springer Line", "L", "List of base lines to create vault springers", GH_ParamAccess.tree);
            pManager.AddBrepParameter("Voussoirs", "V", "Voussoirs to analyse", GH_ParamAccess.tree);
            //pManager.AddPlaneParameter("U Planes", "Pu", "Division Planes with constant U-value.", GH_ParamAccess.tree);
            //pManager.AddNumberParameter("Springer Width", "SpringerWidth", "Distance perpendicular to springer line", GH_ParamAccess.item, 0.3);
        }

        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddBrepParameter("Springers", "S", "Finished Springers", GH_ParamAccess.tree);
            pManager.AddBrepParameter("Voussoirs", "V", "Non transformed voussoirs", GH_ParamAccess.tree);
            //pManager.AddPointParameter("Log", "p1", "All messages generated during execution", GH_ParamAccess.tree);
            //pManager.AddPointParameter("Log", "p2", "All messages generated during execution", GH_ParamAccess.tree);
            //pManager.AddPointParameter("Log", "p3", "All messages generated during execution", GH_ParamAccess.tree);
            //pManager.AddCurveParameter("Voussoirs", "c1", "Non transformed voussoirs", GH_ParamAccess.tree);
            //pManager.AddCurveParameter("Voussoirs", "c2", "Non transformed voussoirs", GH_ParamAccess.tree);
            //pManager.AddCurveParameter("Voussoirs", "c3", "Non transformed voussoirs", GH_ParamAccess.tree);
            //pManager.AddBrepParameter("Springers", "b1", "Finished Springers", GH_ParamAccess.tree);
            //pManager.AddBrepParameter("Springers", "b2", "Finished Springers", GH_ParamAccess.tree);
            //pManager.AddBrepParameter("Springers", "b3", "Finished Springers", GH_ParamAccess.tree);

        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            GH_Structure<GH_Brep> surfacesTree = new GH_Structure<GH_Brep>();
            GH_Structure<GH_Curve> springerLinesTree = new GH_Structure<GH_Curve>();
            GH_Structure<GH_Brep> voussoirsTree = new GH_Structure<GH_Brep>();

            DA.GetDataTree(0, out surfacesTree);
            DA.GetDataTree(1, out springerLinesTree);
            DA.GetDataTree(2, out voussoirsTree);

            //Per Vault
            GH_Structure<GH_Brep> outspringers = new GH_Structure<GH_Brep>();
            GH_Structure<GH_Brep> outvoussoirs = new GH_Structure<GH_Brep>();

            foreach (GH_Path path in surfacesTree.Paths)
            {
                Brep surface = surfacesTree[path][0].Value;
                List<Curve> springerLines = new List<Curve>();
                if (springerLinesTree.PathCount == 1)
                {
                    springerLines = springerLinesTree
                        .AllData(false)
                        .Cast<GH_Curve>()
                        .Select(g => g.Value)
                        .ToList();
                }
                else springerLines = springerLinesTree.Branches[springerLinesTree.Paths.IndexOf(path)].Select(x => x.Value).ToList();

                List<Brep> voussoirs = voussoirsTree.Branches[voussoirsTree.Paths.IndexOf(path)].Select(x => x.Value).ToList();
                var alignedsurface = PolylineUtils.AlignNormalToWorldZ(surface.Surfaces[0]);
                var width = 1;
                // 2️⃣ Compute the centroid of each Brep
                List<Point3d> centroids = voussoirs
                    .Select(v =>
                    {
                        var amp = v.GetBoundingBox(true).Center;
                        return amp;
                    })
                    .ToList();

                List<Brep> restvoussoirs = new List<Brep>();

                // 3️⃣ Compute the average of all centroids → epicenter
                Point3d epicenter = new Point3d(
                    centroids.Average(p => p.X),
                    centroids.Average(p => p.Y),
                    centroids.Average(p => p.Z)
                );

                var line0 = new Line(springerLines[0].PointAtStart, springerLines[0].PointAtEnd);
                var line1 = new Line(springerLines[1].PointAtStart, springerLines[1].PointAtEnd);

                var pl0 = new Plane(line0.From, line0.To, line1.From);

                double maxSize = 0;
                foreach (var b in voussoirs)
                {
                    var bBBox = b.GetBoundingBox(true);
                    var voussSize = bBBox.Diagonal.Length;

                    if (maxSize < voussSize)
                    {
                        maxSize = voussSize;
                    }
                }

                for (int i = 0; i < springerLines.Count; i++)
                {
                    var l = springerLines[i];
                    List<Brep> selectedVoussoirs = new List<Brep>();

                    for (int j = 0; j < voussoirs.Count; j++)
                    {
                        var b = centroids[j];
                        double t;
                        l.ClosestPoint(b, out t);
                        var lp = l.PointAt(t);
                        var dist = lp.DistanceTo(b);

                        if (dist < maxSize / 2)
                        {
                            selectedVoussoirs.Add(voussoirs[j]);
                        }
                    }

                    List<Brep> closestVoussoirs = new List<Brep>();

                    for (int j = 0; j < selectedVoussoirs.Count; j++)
                    {
                        var centroidJ = selectedVoussoirs[j].GetBoundingBox(true).Center;

                        // distance from line to voussoir j
                        double t;
                        l.ClosestPoint(centroidJ, out t);
                        var lp = l.PointAt(t);
                        var C1 = selectedVoussoirs[j].ClosestPoint(lp);
                        var distJ = lp.DistanceTo(C1);

                        bool isclosest = true;

                        for (int k = 0; k < selectedVoussoirs.Count; k++)
                        {
                            if (k == j) continue;
                            var C = selectedVoussoirs[k];
                            var C2 = C.ClosestPoint(lp);
                            var distC = lp.DistanceTo(C2);

                            if (distC < distJ)
                            {
                                isclosest = false;
                                break;
                            }
                        }

                        if (isclosest) closestVoussoirs.Add(selectedVoussoirs[j]);
                    }

                    if (i == 0)
                    {
                        restvoussoirs = voussoirs
                            .Skip(closestVoussoirs.Count)
                            .ToList();
                    }
                    else if (i == 1)
                    {
                        restvoussoirs = restvoussoirs
                            .Take(restvoussoirs.Count - closestVoussoirs.Count)
                            .ToList();
                    }

                    for (int j = 0; j < closestVoussoirs.Count; j++)
                    {
                        var b = closestVoussoirs[j];
                        var bCenter = b.GetBoundingBox(true).Center;

                        double u, v;
                        alignedsurface.ClosestPoint(bCenter, out u, out v);
                        Vector3d srfNormal = alignedsurface.NormalAt(u, v);
                        srfNormal.Unitize();
                        var cPlane = new Plane(bCenter, srfNormal);

                        List<Point3d> vertices = b.Vertices.Select(p => p.Location).ToList();
                        List<Point3d> inVertices = new List<Point3d>();
                        List<Point3d> exVertices = new List<Point3d>();

                        foreach (var p in vertices)
                        {
                            if (cPlane.DistanceTo(p) > 0) exVertices.Add(p);
                            else inVertices.Add(p);
                        }

                        //Debug.WriteLine($"inVertices: {inVertices.Count}; exVertices: {inVertices.Count}");
                        // Sort by Z (lowest → highest)
                        inVertices = inVertices.OrderBy(p => p.Z).ToList();
                        exVertices = exVertices.OrderBy(p => p.Z).ToList();

                        // Split in half
                        int inHalf = inVertices.Count / 2;
                        int exHalf = exVertices.Count / 2;

                        // Lowest half (VL) and highest half (VH)
                        List<Point3d> inVL = inVertices.Take(inHalf).ToList();
                        List<Point3d> inVH = inVertices.Skip(inHalf).ToList();

                        List<Point3d> exVL = exVertices.Take(exHalf).ToList();
                        List<Point3d> exVH = exVertices.Skip(exHalf).ToList();

                        inVL = SprTri01.SortByLine(inVL, l);
                        inVH = SprTri01.SortByLine(inVH, l);
                        exVL = SprTri01.SortByLine(exVL, l);
                        exVH = SprTri01.SortByLine(exVH, l);

                        List<Point3d> face0 = new List<Point3d>();
                        List<Point3d> face1 = new List<Point3d>();

                        face0.Add(inVL[0]);
                        face0.Add(inVH[0]);
                        face0.Add(exVH[0]);
                        face0.Add(exVL[0]);

                        face1.Add(inVL[1]);
                        face1.Add(inVH[1]);
                        face1.Add(exVH[1]);
                        face1.Add(exVL[1]);

                        Plane face0Plane;
                        Plane.FitPlaneToPoints(face0, out face0Plane);
                        Plane face1Plane;
                        Plane.FitPlaneToPoints(face1, out face1Plane);

                        // Face 0
                        // Pt1
                        var e0 = Intersection.CurvePlane(l, face0Plane, 0.1);
                        var pt1f0 = e0[0].PointA;
                        double u1, v1;
                        face0Plane.ClosestParameter(pt1f0, out u1, out v1);
                        pt1f0 = face0Plane.PointAt(u1, v1);

                        // Pt2
                        var pt2f0 = inVH[0];
                        double u2, v2;
                        face0Plane.ClosestParameter(pt2f0, out u2, out v2);
                        pt2f0 = face0Plane.PointAt(u2, v2);

                        // Pt3
                        var pt3f0 = exVH[0];
                        double u21, v21;
                        face0Plane.ClosestParameter(pt3f0, out u21, out v21);
                        pt3f0 = face0Plane.PointAt(u21, v21);

                        // Pt4
                        double z0;
                        var edge0 = new Line(exVH[0], exVL[0]);
                        edge0.Extend(0, maxSize * 10);
                        var e1 = Intersection.LinePlane(edge0, pl0, out z0);
                        var pt4f0 = edge0.PointAt(z0);
                        double u3, v3;
                        face0Plane.ClosestParameter(pt4f0, out u3, out v3);
                        pt4f0 = face0Plane.PointAt(u3, v3);

                        // Face 1
                        // Pt1
                        var e2 = Intersection.CurvePlane(l, face1Plane, 0.1);
                        var pt1f1 = e2[0].PointA;
                        double u4, v4;
                        face1Plane.ClosestParameter(pt1f1, out u4, out v4);
                        pt1f1 = face1Plane.PointAt(u4, v4);

                        // Pt2
                        var pt2f1 = inVH[1];
                        double u5, v5;
                        face1Plane.ClosestParameter(pt2f1, out u5, out v5);
                        pt2f1 = face1Plane.PointAt(u5, v5);

                        // Pt3
                        var pt3f1 = exVH[1];
                        double u51, v51;
                        face1Plane.ClosestParameter(pt3f1, out u51, out v51);
                        pt3f1 = face1Plane.PointAt(u51, v51);

                        // Pt4
                        double z1;
                        var edge1 = new Line(exVH[1], exVL[1]);
                        edge1.Extend(0, maxSize * 10);
                        var e3 = Intersection.LinePlane(edge1, pl0, out z1);
                        var pt4f1 = edge1.PointAt(z1);
                        double u6, v6;
                        face1Plane.ClosestParameter(pt4f1, out u6, out v6);
                        pt4f1 = face1Plane.PointAt(u6, v6);


                        // Build springer faces
                        List<Point3d> face0springer = new List<Point3d>();
                        List<Point3d> face1springer = new List<Point3d>();

                        face0springer.Add(pt1f0);
                        face0springer.Add(pt2f0);
                        face0springer.Add(pt3f0);
                        face0springer.Add(pt4f0);
                        face0springer.Add(pt1f0);

                        face1springer.Add(pt1f1);
                        face1springer.Add(pt2f1);
                        face1springer.Add(pt3f1);
                        face1springer.Add(pt4f1);
                        face1springer.Add(pt1f1);

                        var face0polyline = new Polyline(face0springer);
                        var face1polyline = new Polyline(face1springer);

                        List<Brep> openSpringer = new List<Brep>();
                        for (int p = 0; p < face0springer.Count - 1; p++)
                        {
                            List<Curve> loftlines = new List<Curve>();
                            var lineA = new Line(face0springer[p], face0springer[p + 1]);
                            loftlines.Add(new LineCurve(lineA));
                            var lineB = new Line(face1springer[p], face1springer[p + 1]);
                            loftlines.Add(new LineCurve(lineB));
                            var spface = Brep.CreateFromLoft(loftlines, Point3d.Unset, Point3d.Unset, LoftType.Normal, false);
                            openSpringer.Add(spface[0]);
                        }
                        var springerBody = Brep.JoinBreps(openSpringer, 0.001);

                        var springerface0 = Brep.CreatePlanarBreps(face0polyline.ToNurbsCurve());
                        var springerface1 = Brep.CreatePlanarBreps(face1polyline.ToNurbsCurve());

                        if (springerface0 == null || springerface0.Length == 0 ||
                            springerface1 == null || springerface1.Length == 0)
                        {
                            AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "\nNo springers generated: \nConnection to base already guaranteed.");
                            continue;
                        }

                        var springerbreps = new List<Brep>();
                        springerbreps.Add(springerface0[0]);
                        springerbreps.Add(springerface1[0]);
                        springerbreps.Add(springerBody[0]);

                        var springer = Brep.JoinBreps(springerbreps, 0.1);
                        outspringers.Append(new GH_Brep(springer[0]), path);
                    }

                    //outspringers.AppendRange(closestVoussoirs.Select(p => new GH_Brep(p)));
                }

                outvoussoirs.AppendRange(restvoussoirs.Select(p => new GH_Brep(p)), path);
            }
            DA.SetDataTree(0, outspringers);
            DA.SetDataTree(1, outvoussoirs);
        }
    }
}
