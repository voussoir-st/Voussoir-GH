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
    public static class TreeUtils
    {
        /// <summary>
        /// Partitions each branch of a GH_Structure into sub-branches of a given size.
        /// Equivalent to the Grasshopper "Partition List" component, applied branch-wise.
        /// </summary>
        public static GH_Structure<T> PartitionTree<T>(GH_Structure<T> inputTree, int partitionSize)
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
        public static void CurveClosestPointTree(GH_Structure<GH_Point> pointsTree, GH_Structure<GH_Curve> curvesTree, out GH_Structure<GH_Point> closestPointsTree, out GH_Structure<GH_Number> closestParamsTree,  out GH_Structure<GH_Number> distancesTree)
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
        public static GH_Structure<T> SortListTree<T>(GH_Structure<GH_Number> keysTree, GH_Structure<T> valuesTree) where T : IGH_Goo
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
                //Debug.WriteLine($"branch: " + path);
                if (path.Length > 1)
                {
                    int[] parentIndices = path.Indices.Take(path.Length - 1).ToArray();
                    newPath = new GH_Path(parentIndices);
                }
                else
                {
                    newPath = path;
                }
                //Debug.WriteLine($"newPath: " + newPath);
                foreach (T item in branch.Cast<T>())
                    trimmedTree.Append(item, newPath);
            }

            return trimmedTree;
        }
        public static void SplitTreeAt<T>(GH_Structure<T> inputTree, int splitIndex, out GH_Structure<T> treeA, out GH_Structure<T> treeB) where T : IGH_Goo
        {
            treeA = new GH_Structure<T>();
            treeB = new GH_Structure<T>();

            for (int i = 0; i < inputTree.Branches.Count; i++)
            {
                var path = inputTree.Paths[i];
                var branch = inputTree.Branches[i];

                int count = branch.Count;
                int split = Math.Max(0, Math.Min(splitIndex, count));

                var listA = branch.Take(split).ToList();
                var listB = branch.Skip(split).ToList();

                treeA.AppendRange(listA, path);
                treeB.AppendRange(listB, path);
            }
        }
        public static Brep LoftBySegments(Curve crvA, Curve crvB, double tol = 1e-6)
        {
            // 1. Explode both curves into polyline segments
            var segsA = ExplodeToSegments(crvA);
            var segsB = ExplodeToSegments(crvB);

            //// Ensure same number of segments
            //if (segsA.Count != segsB.Count)
            //    throw new System.Exception("Curves have different segment counts — cannot pair segments.");

            var breps = new List<Brep>();

            // 2. Loft corresponding segments
            for (int i = 0; i < segsA.Count; i++)
            {
                var segA = segsA[i];
                var segB = segsB[i];

                // Create the loft between the two straight segments
                var loft = Brep.CreateFromLoft(
                    new List<Curve> { segA, segB },
                    Point3d.Unset,
                    Point3d.Unset,
                    LoftType.Straight,
                    false
                );

                if (loft != null && loft.Any())
                    breps.Add(loft.First());
            }

            // 3. Optionally join all segment lofts into one Brep
            var joined = Brep.JoinBreps(breps, tol);
            return joined != null && joined.Any() ? joined.First() : null;
        }
        public static List<Curve> ExplodeToSegments(Curve crv)
        {
            // Try to get polyline representation
            if (crv.TryGetPolyline(out Polyline pl))
            {
                var segs = new List<Curve>();
                for (int i = 0; i < pl.Count - 1; i++)
                    segs.Add(new LineCurve(pl[i], pl[i + 1]));
                return segs;
            }

            // Otherwise, explode the curve (handles polylines, polycurves, etc.)
            var exploded = crv.DuplicateSegments()?.ToList();
            if (exploded != null && exploded.Count > 0)
                return exploded;

            // Fallback: divide and create linear approximations
            var pts = crv.DivideEquidistant(10);
            if (pts == null || pts.Length < 2) return new List<Curve>();
            var lines = new List<Curve>();
            for (int i = 0; i < pts.Length - 1; i++)
                lines.Add(new LineCurve(pts[i], pts[i + 1]));
            return lines;
        }


    }
    public class VoussoirCreate : GH_Component
    {
        public VoussoirCreate()
            : base("Springer - Wall", "SprW",
                  "Creates a massive Springer with a flat top",
                  "Voussoir", "Springer") 
        { }

        public override Guid ComponentGuid => new Guid("EC88F9F2-CD3B-4C41-ADFF-FD189794137C");

        protected override System.Drawing.Bitmap Icon => VoussoirPlugin03.Properties.Resources.springer;

        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddBrepParameter("Vault Surface", "VaultSurface", "Surface that defines the vault", GH_ParamAccess.item);
            pManager.AddCurveParameter("Springer Line", "SpringerLine", "List of base lines to create vault springers", GH_ParamAccess.list);
            pManager.AddBrepParameter("Voussoirs", "Voussoirs", "Voussoirs to analyse", GH_ParamAccess.tree);
            pManager.AddPlaneParameter("Transversal Planes", "TransversalPlanes", "Planes at each span division", GH_ParamAccess.list);
            pManager.AddNumberParameter("Springer Width", "SpringerWidth", "Distance perpendicular to springer line", GH_ParamAccess.item, 0.3);
        }

        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddBrepParameter("Springers", "S", "Finished Springers", GH_ParamAccess.tree);
            pManager.AddBrepParameter("Voussoirs", "V", "Non transformed voussoirs", GH_ParamAccess.tree);
            pManager.AddCurveParameter("Top Wall Line", "WL", "Lines to build walls on the springers", GH_ParamAccess.list);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            GH_Structure<GH_Brep> surfacesTree = new GH_Structure<GH_Brep>();
            GH_Structure<GH_Curve> springerLinesTree = new GH_Structure<GH_Curve>();
            GH_Structure<GH_Brep> voussoirsTree = new GH_Structure<GH_Brep>();
            GH_Structure<GH_Plane> transPlanesTree = new GH_Structure<GH_Plane>();

            if (!DA.GetDataTree(0, out surfacesTree)) ;
            if (!DA.GetDataTree(1, out springerLinesTree)) ;
            if (!DA.GetDataTree(2, out voussoirsTree)) ;
            if (!DA.GetDataTree(3, out transPlanesTree)) ;

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

                List<Brep> voussoirs = voussoirsTree.Branches[voussoirsTree.Paths.IndexOf(path)].Select(x => x.Value).ToList();
                List<Plane> transPlanes = transPlanesTree.Branches[transPlanesTree.Paths.IndexOf(path)].Select(x => x.Value).ToList();
                var alignedsurface = PolylineUtils.AlignNormalToWorldZ(surface.Surfaces[0]);
                var line0 = new Line(springerLines[0].PointAtStart, springerLines[0].PointAtEnd);
                var line1 = new Line(springerLines[1].PointAtStart, springerLines[1].PointAtEnd);

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

                double line0t, line1t;
                var splpointAtruth = Intersection.LinePlane(line0, transPlanes[0], out line0t);
                var splpointA = line0.PointAt(line0t);
                var splpointBtruth = Intersection.LinePlane(line1, transPlanes[0], out line1t);
                var splpointB = line1.PointAt(line1t);
                var spllineA = new Line(splpointA, splpointB);
                spllineA.Extend(2, 2);
                c2.Append(new GH_Curve(new LineCurve(spllineA)));

                double line2t, line3t;
                var splpointCtruth = Intersection.LinePlane(line0, transPlanes.Last(), out line2t);
                var splpointC = line0.PointAt(line2t);
                var splpointDtruth = Intersection.LinePlane(line1, transPlanes.Last(), out line3t);
                var splpointD = line1.PointAt(line3t);
                var spllineB = new Line(splpointA, splpointB);
                spllineB.Extend(2, 2);

                //log.Append(new GH_Plane(splinePlane0), path);

                var baseplane0 = new Plane(splpointA, splpointC, spllineA.PointAt(0));
                //pl1.Append(new GH_Plane(baseplane0));
                var baseplane1 = new Plane(splpointB, splpointD, spllineA.PointAt(1));

                Debug.WriteLine($"transPlanes: {transPlanes.Count}");

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
                    List<Brep> intrados = new List<Brep>();
                    List<Brep> extrados = new List<Brep>();

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
                        extrados.Add(faces[1]);
                    }

                    //log.AppendRange(extrados.Select(x => new GH_Brep(x)), path);

                    for (int j = 0; j < springerLines.Count; j++)
                    {

                        if (j == 0)
                        {
                            List<Point3d> springerPolylinePoints1 = new List<Point3d>();
                            List<Point3d> springerPolylinePoints2 = new List<Point3d>();

                            for (int k = 0; k < 2; k++)
                            {
                                var intpointA = Intersection.CurvePlane(springerLines[1], transPlanes[i + k], RhinoMath.ZeroTolerance)[0].PointA;

                                var transformIntrados = intrados[0];
                                //b1.Append(new GH_Brep(transformIntrados));
                                var transformExtrados = extrados[0];
                                //b2.Append(new GH_Brep(transformExtrados));

                                List<Point3d> springerPLinePoints = new List<Point3d>();

                                //Pt 1
                                var pt1 = Intersection.CurvePlane(springerLines[0], transPlanes[i + k], RhinoMath.ZeroTolerance)[0].PointA;
                                //outspringers.Append(new GH_Point(pt1), path);

                                //Pt 2
                                var pt2 = Point3d.Unset;
                                var z = double.MinValue;
                                foreach (var p in transformIntrados.Vertices)
                                {
                                    if (-0.1 < transPlanes[i + k].DistanceTo(p.Location) && transPlanes[i + k].DistanceTo(p.Location) < 0.1)
                                    {
                                        if (p.Location.Z > z)
                                        {
                                            pt2 = p.Location;
                                            z = p.Location.Z;
                                        }
                                    }
                                }
                                //outspringers.Append(new GH_Point(pt2), path);

                                //Pt 3
                                var pt3 = Point3d.Unset;
                                var ez = double.MinValue;
                                foreach (var p in transformExtrados.Vertices)
                                {
                                    if (-0.01 < transPlanes[i + k].DistanceTo(p.Location) && transPlanes[i + k].DistanceTo(p.Location) < 0.01)
                                    {
                                        if (p.Location.Z > ez)
                                        {
                                            pt3 = p.Location;
                                            ez = p.Location.Z;
                                        }
                                    }
                                }
                                //outspringers.Append(new GH_Point(pt3), path);
                                List<BrepVertex> extradospoints = extrados[0].Vertices.ToList();
                                //outspringers.AppendRange(extradospoints.Select(x => new GH_Point(x.Location)), path);
                                //p1.Append(new GH_Point(pt3));
                                //Pt 4                                
                                var lpoint = Point3d.Unset;
                                var lz = double.MaxValue;

                                foreach (var p in transformExtrados.Vertices)
                                {
                                    //p3.Append(new GH_Point(p.Location));
                                    if (-0.1 < transPlanes[i + k].DistanceTo(p.Location) && transPlanes[i + k].DistanceTo(p.Location) < 0.1)
                                    {
                                        if (p.Location.Z < lz)
                                        {
                                            lpoint = p.Location;
                                            lz = p.Location.Z;
                                        }
                                    }
                                }

                                var dirLine = new Line(pt3, lpoint);
                                dirLine.Extend(300, 300);
                                var intLine = new Line(pt1, intpointA);
                                intLine.Extend(300, 300);

                                var pt4 = Point3d.Unset;
                                double t;
                                var truth = Intersection.LinePlane(dirLine, baseplane0, out t);
                                pt4 = dirLine.PointAt(t);

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
                            Polyline pl1 = new Polyline(springerPolylinePoints1);
                            Polyline pl2 = new Polyline(springerPolylinePoints2);
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
                                var intpointA = Intersection.CurvePlane(springerLines[0], transPlanes[i + k], RhinoMath.ZeroTolerance)[0].PointA;

                                var transformIntrados = intrados.Last();
                                //b1.Append(new GH_Brep(transformIntrados));
                                var transformExtrados = extrados.Last();
                                //b2.Append(new GH_Brep(transformExtrados));

                                List<Point3d> springerPLinePoints = new List<Point3d>();

                                //Pt 1
                                var pt1 = Intersection.CurvePlane(springerLines[1], transPlanes[i + k], RhinoMath.ZeroTolerance)[0].PointA;
                                //outspringers.Append(new GH_Point(pt1), path);

                                //Pt 2
                                var pt2 = Point3d.Unset;
                                var z = double.MinValue;
                                foreach (var p in transformIntrados.Vertices)
                                {
                                    if (-0.1 < transPlanes[i + k].DistanceTo(p.Location) && transPlanes[i + k].DistanceTo(p.Location) < 0.1)
                                    {
                                        if (p.Location.Z > z)
                                        {
                                            pt2 = p.Location;
                                            z = p.Location.Z;
                                        }
                                    }
                                }
                                //outspringers.Append(new GH_Point(pt2), path);

                                //Pt 3
                                var pt3 = Point3d.Unset;
                                var ez = double.MinValue;
                                foreach (var p in transformExtrados.Vertices)
                                {
                                    if (-0.01 < transPlanes[i + k].DistanceTo(p.Location) && transPlanes[i + k].DistanceTo(p.Location) < 0.01)
                                    {
                                        if (p.Location.Z > ez)
                                        {
                                            pt3 = p.Location;
                                            ez = p.Location.Z;
                                        }
                                    }
                                }
                                //outspringers.Append(new GH_Point(pt3), path);
                                List<BrepVertex> extradospoints = extrados[0].Vertices.ToList();
                                //outspringers.AppendRange(extradospoints.Select(x => new GH_Point(x.Location)), path);
                                //p1.Append(new GH_Point(pt3));
                                //Pt 4                                
                                var lpoint = Point3d.Unset;
                                var lz = double.MaxValue;

                                foreach (var p in transformExtrados.Vertices)
                                {
                                    //p3.Append(new GH_Point(p.Location));
                                    if (-0.1 < transPlanes[i + k].DistanceTo(p.Location) && transPlanes[i + k].DistanceTo(p.Location) < 0.1)
                                    {
                                        if (p.Location.Z < lz)
                                        {
                                            lpoint = p.Location;
                                            lz = p.Location.Z;
                                        }
                                    }
                                }

                                var dirLine = new Line(pt3, lpoint);
                                dirLine.Extend(300, 300);
                                var intLine = new Line(pt1, intpointA);
                                intLine.Extend(300, 300);

                                var pt4 = Point3d.Unset;
                                double t;
                                var truth = Intersection.LinePlane(dirLine, baseplane1, out t);
                                pt4 = dirLine.PointAt(t);

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
                            Polyline pl1 = new Polyline(springerPolylinePoints1);
                            Polyline pl2 = new Polyline(springerPolylinePoints2);
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
