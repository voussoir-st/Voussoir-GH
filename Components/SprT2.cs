using Components;
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
using System.Windows.Forms;
using static Rhino.DocObjects.PhysicallyBasedMaterial;

namespace VoussoirPlugin03.Components
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
            pManager.AddBrepParameter("Vault Surface", "VaultSurface", "Surface that defines the vault", GH_ParamAccess.tree);
            pManager.AddCurveParameter("Springer Line", "SpringerLine", "List of base lines to create vault springers", GH_ParamAccess.tree);
            pManager.AddBrepParameter("Voussoirs", "Voussoirs", "Voussoirs to analyse", GH_ParamAccess.tree);
            pManager.AddPlaneParameter("Transversal Planes", "TransversalPlanes", "Planes at each span division", GH_ParamAccess.tree);
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
            GH_Structure<GH_Plane> transPlanesTree = new GH_Structure<GH_Plane>();

            if (!DA.GetDataTree(0, out surfacesTree));
            if (!DA.GetDataTree(1, out springerLinesTree));
            if (!DA.GetDataTree(2, out voussoirsTree));
            if (!DA.GetDataTree(3, out transPlanesTree));

            //Per Vault
            GH_Structure<GH_Brep> outspringers = new GH_Structure<GH_Brep>();
            GH_Structure<GH_Brep> outvoussoirs = new GH_Structure<GH_Brep>();
            GH_Structure<GH_Point> p1 = new GH_Structure<GH_Point>();
            GH_Structure<GH_Point> p2 = new GH_Structure<GH_Point>();
            GH_Structure<GH_Point> p3 = new GH_Structure<GH_Point>();
            GH_Structure<GH_Curve> c1 = new GH_Structure<GH_Curve>();
            GH_Structure<GH_Curve> c2 = new GH_Structure<GH_Curve>();
            GH_Structure<GH_Curve> c3 = new GH_Structure<GH_Curve>();
            GH_Structure<GH_Brep> b1 = new GH_Structure<GH_Brep>();
            GH_Structure<GH_Brep> b2 = new GH_Structure<GH_Brep>();
            GH_Structure<GH_Brep> b3 = new GH_Structure<GH_Brep>();

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
                var width = 1;
                List<Brep> voussoirs = voussoirsTree.Branches[voussoirsTree.Paths.IndexOf(path)].Select(x => x.Value).ToList();
                List<Plane> transPlanes = transPlanesTree.Branches[transPlanesTree.Paths.IndexOf(path)].Select(x => x.Value).ToList();
                var alignedsurface = PolylineUtils.AlignNormalToWorldZ(surface.Surfaces[0]);
               
                // 2️⃣ Compute the centroid of each Brep
                List<Point3d> centroids = voussoirs
                    .Select(v =>
                    {
                        // Use area mass properties to get centroid
                        var amp = AreaMassProperties.Compute(v);
                        return amp.Centroid;
                    })
                    .ToList();

                // 3️⃣ Compute the average of all centroids → epicenter
                Point3d epicenter = new Point3d(
                    centroids.Average(p => p.X),
                    centroids.Average(p => p.Y),
                    centroids.Average(p => p.Z)
                );

                var line0 = new Line(springerLines[0].PointAtStart, springerLines[0].PointAtEnd);
                var line1 = new Line(springerLines[1].PointAtStart, springerLines[1].PointAtEnd);

                var pl0 = new Plane(line0.From, line0.To, line0.From + Plane.WorldXY.ZAxis * 1);
                var zdline = new Line(epicenter, line0.PointAt(0.5)).Direction;
                zdline.Unitize();
                if (pl0.ZAxis * zdline < 0)
                    pl0.Flip();

                var line2 = new Line(line0.From + pl0.ZAxis * width, line0.To + pl0.ZAxis * width);
                line2.Extend(300, 300);

                var pla1 = new Plane(line1.From, line1.To, line1.From + Plane.WorldXY.ZAxis * 1);
                var zdline1 = new Line(epicenter, line1.PointAt(0.5)).Direction;
                zdline1.Unitize();
                if (pla1.ZAxis * zdline < 0)
                    pla1.Flip();

                var line3 = new Line(line0.From + pl0.ZAxis * width, line0.To + pl0.ZAxis * width);
                line3.Extend(300, 300);

                List<Vector3d> intLines = new List<Vector3d>();

                foreach (var p in transPlanes)
                {
                    double a, b;
                    Intersection.LinePlane(line0, p, out a);
                    var c = line0.PointAt(a);
                    Intersection.LinePlane(line2, p, out b);
                    var d = line2.PointAt(b);
                    var e = new Line(c, d).Direction;
                    e.Unitize();
                    intLines.Add(e);
                }
                var baseplane0 = new Plane(line0.From, line0.To, line0.From + intLines[0] * 1);
                var baseplane1 = new Plane(line1.From, line1.To, line1.From + -intLines[0] * 1);

                Debug.WriteLine($"transPlanes: {transPlanes.Count}");
                Debug.WriteLine($"transPlanes: {transPlanes.Count}");
                List<Point3d> pts40 = new List<Point3d>();
                List<Point3d> pts41 = new List<Point3d>();
                List<Brep> intrados = new List<Brep>();
                List<Brep> extrados = new List<Brep>();

                //Per Row
                for (int i = 0; i < transPlanes.Count - 1; i++)
                {
                    var intplaneoriginA = transPlanes[i].Origin;
                    var intplaneoriginB = transPlanes[i + 1].Origin;
                    var avgorigin = new Point3d(
                        (intplaneoriginA.X + intplaneoriginB.X) / 2,
                        (intplaneoriginA.Y + intplaneoriginB.Y) / 2,
                        (intplaneoriginA.Z + intplaneoriginB.Z) / 2
                        );

                    var avgPlane = new Plane(avgorigin, transPlanes[i].XAxis, transPlanes[i].YAxis);

                    var planesTol = intplaneoriginA.DistanceTo(intplaneoriginB) * 0.8;

                    List<Brep> rowVoussoirs = new List<Brep>();

                    foreach (Brep v in voussoirs)
                    {
                        Point3d avgVoussoirPt = new Point3d(
                            v.Vertices.Average(p => p.Location.X),
                            v.Vertices.Average(p => p.Location.Y),
                            v.Vertices.Average(p => p.Location.Z)
                        );

                        if (-planesTol / 2 < avgPlane.DistanceTo(avgVoussoirPt) && avgPlane.DistanceTo(avgVoussoirPt) < planesTol / 2)
                        {
                            rowVoussoirs.Add(v);
                        }
                    }
                    Debug.WriteLine($"rowVoussoirs: {rowVoussoirs.Count}");

                    Point3d pointA = Intersection.CurvePlane(springerLines[0], transPlanes[i], RhinoMath.ZeroTolerance)[0].PointA;
                    Point3d pointB = Intersection.CurvePlane(springerLines[1], transPlanes[i], RhinoMath.ZeroTolerance)[0].PointA;

                    var orderLine = new Line(pointA, pointB);

                    var orderedRowVoussoirs = rowVoussoirs.OrderBy(v =>
                    {
                        double t;

                        var pts = v.Vertices.Select(bv => bv.Location).ToList();
                        Point3d p = new Point3d(
                            pts.Average(pt => pt.X),
                            pts.Average(pt => pt.Y),
                            pts.Average(pt => pt.Z)
                        );

                        t = orderLine.ClosestParameter(p);
                        return t;
                    }).ToList();
                    List<Brep> culledvoussoirs = new List<Brep>(orderedRowVoussoirs);
                    culledvoussoirs.RemoveAt(orderedRowVoussoirs.Count - 1);
                    culledvoussoirs.RemoveAt(0);
                    outvoussoirs.AppendRange(culledvoussoirs.Select(v => new GH_Brep(v)), path);

                    //Intrados - Extrados
                    List<Brep> rowintrados = new List<Brep>();
                    List<Brep> rowextrados = new List<Brep>();
                    foreach (Brep v in orderedRowVoussoirs)
                    {
                        //outspringers.Append(new GH_Brep(v), path);
                        // Duplicate the faces of the voussoir
                        List<Brep> faces = v.Faces.Select(f => f.DuplicateFace(true)).ToList();

                        // Get vertices and compute the centroid
                        List<Point3d> points = v.Vertices.Select(p => p.Location).ToList();
                        Point3d voussCenter = new Point3d(
                            points.Average(pt => pt.X),
                            points.Average(pt => pt.Y),
                            points.Average(pt => pt.Z)
                        );

                        // Normal of the aligned surface at centroid
                        Vector3d surfaceNormal = PolylineUtils.GetBrepNormalAtPoint(alignedsurface.ToBrep(), voussCenter, RhinoMath.ZeroTolerance);
                        surfaceNormal.Unitize();

                        // Compute dot products for each face
                        List<double> products = faces
                            .Select(f =>
                            {
                                var pointsf = f.Vertices.Select(vx => vx.Location).ToList();
                                Point3d faceCenter = new Point3d(
                                    points.Average(pt => pt.X),
                                    points.Average(pt => pt.Y),
                                    points.Average(pt => pt.Z)
                                );
                                // Normal of the face at centroid
                                Vector3d faceNormal = PolylineUtils.GetBrepNormalAtPoint(f, faceCenter, RhinoMath.ZeroTolerance);
                                faceNormal.Unitize();

                                // Dot product
                                return faceNormal * surfaceNormal;
                            })
                            .ToList();

                        // Reorder faces according to products (ascending)
                        faces = faces
                            .Select((f, l) => new { Face = f, Product = products[l] })
                            .OrderBy(x => x.Product)
                            .Select(x => x.Face)
                            .ToList();

                        intrados.Add(faces[0]);
                        rowintrados.Add(faces[0]);
                        extrados.Add(faces[1]);
                        rowextrados.Add(faces[1]);

                    }

                    for (int j = 0; j < springerLines.Count; j++)
                    {
                        if (j == 0)
                        {
                            for (int k = 0; k < 2; k++)
                            {
                                //Pt 1
                                var pt1 = Intersection.CurvePlane(springerLines[0], transPlanes[i + k], RhinoMath.ZeroTolerance)[0].PointA;

                                var intpointA = Intersection.CurvePlane(springerLines[1], transPlanes[i + k], RhinoMath.ZeroTolerance)[0].PointA;

                                var plane = transPlanes[i + k];
                                var transformIntrados = rowintrados[0];
                                var transformExtrados = rowextrados[0];
                                //Pt 3
                                var pt3 = transformExtrados.Vertices
                                        .Select(v => v.Location)
                                        .OrderBy(p => Math.Abs(plane.DistanceTo(p))) // closeness to plane
                                        .Take(2)                                     // two closest
                                        .OrderByDescending(p => p.Z)                 // highest Z
                                        .First();
                                //p2.Append(new GH_Point(pt3));

                                //Pt 4                                
                                var lpoint = transformExtrados.Vertices
                                    .Select(v => v.Location)
                                    .OrderBy(p => Math.Abs(plane.DistanceTo(p)))
                                    .Take(2)
                                    .OrderBy(p => p.Z)
                                    .First();
                                //p3.Append(new GH_Point(lpoint));

                                var dirLine = new Line(pt3, lpoint);
                                dirLine.Extend(300, 300);
                                var intLine = new Line(pt1, intpointA);
                                intLine.Extend(300, 300);

                                var pt4 = Point3d.Unset;
                                double t;
                                var truth = Intersection.LinePlane(dirLine, baseplane0, out t);
                                pt4 = dirLine.PointAt(t);
                                //p1.Append(new GH_Point(pt4));
                                pts40.Add(pt4);
                            }
                        }

                        else
                        {
                            for (int k = 0; k < 2; k++)
                            {
                                //Pt 1
                                var pt1 = Intersection.CurvePlane(springerLines[0], transPlanes[i + k], RhinoMath.ZeroTolerance)[0].PointA;

                                var intpointA = Intersection.CurvePlane(springerLines[1], transPlanes[i + k], RhinoMath.ZeroTolerance)[0].PointA;

                                var plane = transPlanes[i + k];
                                var transformIntrados = rowintrados.Last();
                                var transformExtrados = rowextrados.Last();
                                //Pt 3
                                var pt3 = transformExtrados.Vertices
                                        .Select(v => v.Location)
                                        .OrderBy(p => Math.Abs(plane.DistanceTo(p))) // closeness to plane
                                        .Take(2)                                     // two closest
                                        .OrderByDescending(p => p.Z)                 // highest Z
                                        .First();
                                //p2.Append(new GH_Point(pt3));

                                //Pt 4                                
                                var lpoint = transformExtrados.Vertices
                                    .Select(v => v.Location)
                                    .OrderBy(p => Math.Abs(plane.DistanceTo(p)))
                                    .Take(2)
                                    .OrderBy(p => p.Z)
                                    .First();
                                //p3.Append(new GH_Point(lpoint));

                                var dirLine = new Line(pt3, lpoint);
                                dirLine.Extend(300, 300);
                                var intLine = new Line(pt1, intpointA);
                                intLine.Extend(300, 300);



                                var pt4 = Point3d.Unset;
                                double t;
                                var truth = Intersection.LinePlane(dirLine, baseplane1, out t);
                                pt4 = dirLine.PointAt(t);
                                pts41.Add(pt4);
                            }
                        }
                    }
                }

                var Line20 = new Line(pts40[0], pts40.Last());
                c1.Append(new GH_Curve(new LineCurve(Line20)));
                var Line21 = new Line(pts41[0], pts41.Last());
                c1.Append(new GH_Curve(new LineCurve(Line21)));

                for (int i = 0; i < transPlanes.Count - 1; i++)
                {
                    var intplaneoriginA = transPlanes[i].Origin;
                    var intplaneoriginB = transPlanes[i + 1].Origin;
                    var avgorigin = new Point3d(
                        (intplaneoriginA.X + intplaneoriginB.X) / 2,
                        (intplaneoriginA.Y + intplaneoriginB.Y) / 2,
                        (intplaneoriginA.Z + intplaneoriginB.Z) / 2
                        );

                    var avgPlane = new Plane(avgorigin, transPlanes[i].XAxis, transPlanes[i].YAxis);

                    var planesTol = intplaneoriginA.DistanceTo(intplaneoriginB) * 0.8;

                    List<Brep> rowIntrados = new List<Brep>();
                    List<Brep> rowExtrados = new List<Brep>();

                    foreach (Brep v in intrados)
                    {
                        Point3d avgintradosPt = new Point3d(
                            v.Vertices.Average(p => p.Location.X),
                            v.Vertices.Average(p => p.Location.Y),
                            v.Vertices.Average(p => p.Location.Z)
                        );

                        if (-planesTol / 2 < avgPlane.DistanceTo(avgintradosPt) && avgPlane.DistanceTo(avgintradosPt) < planesTol / 2)
                        {
                            rowIntrados.Add(v);
                        }
                    }
                    foreach (Brep v in extrados)
                    {
                        Point3d avgextradosPt = new Point3d(
                            v.Vertices.Average(p => p.Location.X),
                            v.Vertices.Average(p => p.Location.Y),
                            v.Vertices.Average(p => p.Location.Z)
                        );

                        if (-planesTol / 2 < avgPlane.DistanceTo(avgextradosPt) && avgPlane.DistanceTo(avgextradosPt) < planesTol / 2)
                        {
                            rowExtrados.Add(v);
                        }
                    }

                    for (int j = 0; j < springerLines.Count; j++)
                    {

                        if (j == 0)
                        {
                            List<Point3d> springerPolylinePoints1 = new List<Point3d>();
                            List<Point3d> springerPolylinePoints2 = new List<Point3d>();

                            for (int k = 0; k < 2; k++)
                            {
                                var plane = transPlanes[i + k];
                                //pl1.Append(new GH_Plane(plane));
                                var intpointA = Intersection.CurvePlane(springerLines[1], transPlanes[i + k], RhinoMath.ZeroTolerance)[0].PointA;

                                var transformIntrados = rowIntrados[0];
                                //b1.Append(new GH_Brep(transformIntrados));
                                var transformExtrados = rowExtrados[0];
                                //b2.Append(new GH_Brep(transformExtrados));

                                List<Point3d> springerPLinePoints = new List<Point3d>();

                                //Pt 1
                                var pt1 = Intersection.CurvePlane(springerLines[0], transPlanes[i + k], RhinoMath.ZeroTolerance)[0].PointA;
                                p1.Append(new GH_Point(pt1));
                                //outspringers.Append(new GH_Point(pt1), path);
                                //p1.Append(new GH_Point(pt1));

                                //Pt 2
                                var pt2 = transformIntrados.Vertices
                                    .Select(v => v.Location)
                                    .OrderBy(p => Math.Abs(plane.DistanceTo(p))) // closeness to plane
                                    .Take(2)                                     // two closest
                                    .OrderByDescending(p => p.Z)                 // highest Z
                                    .First();

                                //Pt 3
                                var pt3 = transformExtrados.Vertices
                                    .Select(v => v.Location)
                                    .OrderBy(p => Math.Abs(plane.DistanceTo(p))) // closeness to plane
                                    .Take(2)                                     // two closest
                                    .OrderByDescending(p => p.Z)                 // highest Z
                                    .First();
                                //p2.Append(new GH_Point(pt3));

                                //Pt 4                                
                                var lpoint = transformExtrados.Vertices
                                    .Select(v => v.Location)
                                    .OrderBy(p => Math.Abs(plane.DistanceTo(p)))
                                    .Take(2)
                                    .OrderBy(p => p.Z)
                                    .First();
                                //p3.Append(new GH_Point(lpoint));
                                p2.Append(new GH_Point(lpoint));
                                var dirLine = new Line(pt3, lpoint);
                                dirLine.Extend(300, 300);

                                var pt4 = Point3d.Unset;
                                double t;
                                Intersection.LinePlane(Line20, plane, out t);
                                pt4 = Line20.PointAt(t);

                                p3.Append(new GH_Point(pt4));

                                //c1.Append(new GH_Curve(new LineCurve(dirLine)), path);
                                //outspringers.Append(new GH_Point(pt4), path);
                                if (k == 0)
                                {
                                    springerPolylinePoints1.Add(pt1);
                                    springerPolylinePoints1.Add(pt2);
                                    springerPolylinePoints1.Add(pt3);
                                    springerPolylinePoints1.Add(pt4);
                                    springerPolylinePoints1.Add(pt1);
                                }
                                else
                                {
                                    springerPolylinePoints2.Add(pt1);
                                    springerPolylinePoints2.Add(pt2);
                                    springerPolylinePoints2.Add(pt3);
                                    springerPolylinePoints2.Add(pt4);
                                    springerPolylinePoints2.Add(pt1);
                                }
                            }

                            // Assume springerPolylinePoints1 and springerPolylinePoints2 are your point lists
                            Polyline pol1 = new Polyline(springerPolylinePoints1);
                            Polyline pol2 = new Polyline(springerPolylinePoints2);
                            List<Brep> openSpringer = new List<Brep>();
                            for (int p = 0; p < springerPolylinePoints1.Count - 1; p++)
                            {
                                List<Curve> loftlines = new List<Curve>();
                                var lineA = new Line(springerPolylinePoints1[p], springerPolylinePoints1[p + 1]);
                                loftlines.Add(new LineCurve(lineA));
                                var lineB = new Line(springerPolylinePoints2[p], springerPolylinePoints2[p + 1]);
                                loftlines.Add(new LineCurve(lineB));
                                var spface = Brep.CreateFromLoft(loftlines, Point3d.Unset, Point3d.Unset, LoftType.Normal, false);
                                openSpringer.Add(spface[0]);
                            }
                            var springerBody = Brep.JoinBreps(openSpringer, 0.001);

                            // Create planar caps for top and bottom
                            var face1 = NurbsSurface.CreateFromCorners(springerPolylinePoints1[0], springerPolylinePoints1[1], springerPolylinePoints1[2], springerPolylinePoints1[3]);
                            var face2 = NurbsSurface.CreateFromCorners(springerPolylinePoints2[0], springerPolylinePoints2[1], springerPolylinePoints2[2], springerPolylinePoints2[3]);

                            // Join the loft body with the planar caps
                            List<Brep> springerParts = new List<Brep>();
                            springerParts.Add(face1.ToBrep());
                            springerParts.Add(springerBody[0]);
                            springerParts.Add(face2.ToBrep());

                            Brep[] joinedSpringer = Brep.JoinBreps(springerParts, 0.001);

                            outspringers.Append(new GH_Brep(joinedSpringer[0]), path);
                            //log.AppendRange(
                            //    springerPolylinePoints1.Select(pt => new GH_Point(pt)).ToList(),
                            //    path
                            //);
                        }

                        else
                        {
                            List<Point3d> springerPolylinePoints1 = new List<Point3d>();
                            List<Point3d> springerPolylinePoints2 = new List<Point3d>();

                            for (int k = 0; k < 2; k++)
                            {
                                var plane = transPlanes[i + k];

                                var intpointA = Intersection.CurvePlane(springerLines[0], transPlanes[i + k], RhinoMath.ZeroTolerance)[0].PointA;

                                var transformIntrados = rowIntrados.Last();
                                //b1.Append(new GH_Brep(transformIntrados));
                                var transformExtrados = rowExtrados.Last();
                                //b2.Append(new GH_Brep(transformExtrados));

                                List<Point3d> springerPLinePoints = new List<Point3d>();

                                //Pt 1
                                var pt1 = Intersection.CurvePlane(springerLines[1], transPlanes[i + k], RhinoMath.ZeroTolerance)[0].PointA;
                                //outspringers.Append(new GH_Point(pt1), path);

                                //Pt 3
                                var pt2 = transformIntrados.Vertices
                                    .Select(v => v.Location)
                                    .OrderBy(p => Math.Abs(plane.DistanceTo(p))) // closeness to plane
                                    .Take(2)                                     // two closest
                                    .OrderByDescending(p => p.Z)                 // highest Z
                                    .First();
                                //Pt 3
                                var pt3 = transformExtrados.Vertices
                                    .Select(v => v.Location)
                                    .OrderBy(p => Math.Abs(plane.DistanceTo(p))) // closeness to plane
                                    .Take(2)                                     // two closest
                                    .OrderByDescending(p => p.Z)                 // highest Z
                                    .First();

                                //Pt 4                                
                                var lpoint = transformExtrados.Vertices
                                    .Select(v => v.Location)
                                    .OrderBy(p => Math.Abs(plane.DistanceTo(p)))
                                    .Take(2)
                                    .OrderBy(p => p.Z)
                                    .First();

                                var dirLine = new Line(pt3, lpoint);
                                dirLine.Extend(300, 300);
                                var intLine = new Line(pt1, intpointA);
                                intLine.Extend(300, 300);

                                var pt4 = Point3d.Unset;
                                double t;
                                Intersection.LinePlane(Line21, plane, out t);
                                pt4 = Line21.PointAt(t);

                                //p2.Append(new GH_Point(lpoint));

                                //c1.Append(new GH_Curve(new LineCurve(dirLine)), path);
                                //outspringers.Append(new GH_Point(pt4), path);
                                if (k == 0)
                                {
                                    springerPolylinePoints1.Add(pt1);
                                    springerPolylinePoints1.Add(pt2);
                                    springerPolylinePoints1.Add(pt3);
                                    springerPolylinePoints1.Add(pt4);
                                    springerPolylinePoints1.Add(pt1);
                                }
                                else
                                {
                                    springerPolylinePoints2.Add(pt1);
                                    springerPolylinePoints2.Add(pt2);
                                    springerPolylinePoints2.Add(pt3);
                                    springerPolylinePoints2.Add(pt4);
                                    springerPolylinePoints2.Add(pt1);
                                }
                            }

                            // Assume springerPolylinePoints1 and springerPolylinePoints2 are your point lists
                            Polyline pol1 = new Polyline(springerPolylinePoints1);
                            Polyline pol2 = new Polyline(springerPolylinePoints2);
                            List<Brep> openSpringer = new List<Brep>();
                            for (int p = 0; p < springerPolylinePoints1.Count - 1; p++)
                            {
                                List<Curve> loftlines = new List<Curve>();
                                var lineA = new Line(springerPolylinePoints1[p], springerPolylinePoints1[p + 1]);
                                loftlines.Add(new LineCurve(lineA));
                                var lineB = new Line(springerPolylinePoints2[p], springerPolylinePoints2[p + 1]);
                                loftlines.Add(new LineCurve(lineB));
                                var spface = Brep.CreateFromLoft(loftlines, Point3d.Unset, Point3d.Unset, LoftType.Normal, false);
                                openSpringer.Add(spface[0]);
                            }
                            var springerBody = Brep.JoinBreps(openSpringer, 0.001);

                            // Create planar caps for top and bottom
                            var face1 = NurbsSurface.CreateFromCorners(springerPolylinePoints1[0], springerPolylinePoints1[1], springerPolylinePoints1[2], springerPolylinePoints1[3]);
                            var face2 = NurbsSurface.CreateFromCorners(springerPolylinePoints2[0], springerPolylinePoints2[1], springerPolylinePoints2[2], springerPolylinePoints2[3]);

                            // Join the loft body with the planar caps
                            List<Brep> springerParts = new List<Brep>();
                            springerParts.Add(face1.ToBrep());
                            springerParts.Add(springerBody[0]);
                            springerParts.Add(face2.ToBrep());

                            Brep[] joinedSpringer = Brep.JoinBreps(springerParts, 0.001);

                            outspringers.Append(new GH_Brep(joinedSpringer[0]), path);
                        }
                    }
                }
            }
            DA.SetDataTree(0, outspringers);
            DA.SetDataTree(1, outvoussoirs);
            //DA.SetDataTree(2, p1);
            //DA.SetDataTree(3, p2);
            //DA.SetDataTree(4, p3);
            //DA.SetDataTree(5, c1);
            //DA.SetDataTree(6, c2);
            //DA.SetDataTree(7, c3);
            //DA.SetDataTree(8, b1);
            //DA.SetDataTree(9, b2);
            //DA.SetDataTree(10, b3);

        }
    }
}
