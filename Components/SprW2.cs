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
            pManager.AddBrepParameter("Vault Surface", "VaultSurface", "Surface that defines the vault", GH_ParamAccess.tree);
            pManager.AddCurveParameter("Springer Line", "SpringerLine", "List of base lines to create vault springers", GH_ParamAccess.tree);
            pManager.AddBrepParameter("Voussoirs", "Voussoirs", "Voussoirs to analyse", GH_ParamAccess.tree);
            pManager.AddPlaneParameter("Transversal Planes", "TransversalPlanes", "Planes at each span division", GH_ParamAccess.tree);
            pManager.AddNumberParameter("Springer Width", "SpringerWidth", "Distance perpendicular to springer line", GH_ParamAccess.tree, 0.3);
        }

        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddBrepParameter("Springers", "S", "Finished Springers", GH_ParamAccess.tree);
            pManager.AddBrepParameter("Voussoirs", "V", "Non transformed voussoirs", GH_ParamAccess.tree);
            pManager.AddCurveParameter("Top Wall Line", "WL", "Lines to build walls on the springers", GH_ParamAccess.tree);
        //    pManager.AddPointParameter("Log", "p1", "All messages generated during execution", GH_ParamAccess.tree);
        //    pManager.AddPointParameter("Log", "p2", "All messages generated during execution", GH_ParamAccess.tree);
        //    pManager.AddPointParameter("Log", "p3", "All messages generated during execution", GH_ParamAccess.tree);
        //    pManager.AddCurveParameter("Voussoirs", "c1", "Non transformed voussoirs", GH_ParamAccess.tree);
        //    pManager.AddCurveParameter("Voussoirs", "c2", "Non transformed voussoirs", GH_ParamAccess.tree);
        //    pManager.AddCurveParameter("Voussoirs", "c3", "Non transformed voussoirs", GH_ParamAccess.tree);
        //    pManager.AddBrepParameter("Springers", "b1", "Finished Springers", GH_ParamAccess.tree);
        //    pManager.AddBrepParameter("Springers", "b2", "Finished Springers", GH_ParamAccess.tree);
        //    pManager.AddBrepParameter("Springers", "b3", "Finished Springers", GH_ParamAccess.tree);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            GH_Structure<GH_Brep> surfacesTree = new GH_Structure<GH_Brep>();
            GH_Structure<GH_Curve> springerLinesTree = new GH_Structure<GH_Curve>();
            GH_Structure<GH_Brep> voussoirsTree = new GH_Structure<GH_Brep>();
            GH_Structure<GH_Plane> transPlanesTree = new GH_Structure<GH_Plane>();
            GH_Structure<GH_Number> springerWidthTree = new GH_Structure<GH_Number>();

            if (!DA.GetDataTree(0, out surfacesTree));
            if (!DA.GetDataTree(1, out springerLinesTree));
            if (!DA.GetDataTree(2, out voussoirsTree));
            if (!DA.GetDataTree(3, out transPlanesTree));
            if (!DA.GetDataTree(4, out springerWidthTree)) ;

            //Per Vault
            GH_Structure<GH_Brep> outspringers = new GH_Structure<GH_Brep>();
            GH_Structure<GH_Brep> outvoussoirs = new GH_Structure<GH_Brep>();
            GH_Structure<GH_Curve> WL = new GH_Structure<GH_Curve>();

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

                double width = 0.3;
                if (springerWidthTree.Branches.Count == 1)
                    width = springerWidthTree[0][0].Value;
                else width = springerWidthTree[path][0].Value;

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

                var pl1 = new Plane(line1.From, line1.To, line1.From + Plane.WorldXY.ZAxis * 1);
                var zdline1 = new Line(epicenter, line1.PointAt(0.5)).Direction;
                zdline1.Unitize();
                if (pl1.ZAxis * zdline < 0)
                    pl1.Flip();

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

                Debug.WriteLine($"transPlanes: {transPlanes.Count}");

                List<int> side0 = new List<int>();
                List<int> side1 = new List<int>();
                List<Brep> side0Voussoirs = new List<Brep>();
                List<Brep> side1Voussoirs = new List<Brep>();


                List<Point3d> p1 = new List<Point3d>();
                List<Point3d> p2 = new List<Point3d>();

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

                    List<Brep> selectedvoussoirs1 = new List<Brep>();
                    List<Brep> selectedvoussoirs2 = new List<Brep>();

                    for (int j = 0; j < 2; j++)
                    {
                        if (j == 0)
                        {
                            var vertplane = new Plane(springerLines[0].PointAtStart, springerLines[0].PointAtEnd, springerLines[0].PointAtStart + Plane.WorldXY.ZAxis * 1);
                            var prez = vertplane.ZAxis;
                            prez.Unitize();

                            var dLine = new Line(epicenter, line0.PointAt(0.5)).Direction;
                            dLine.Unitize();
                            if (prez * dLine < 0) vertplane.Flip();

                            foreach (Brep v in orderedRowVoussoirs)
                            {
                                List<Point3d> voussoir = new List<Point3d>();
                                foreach (var p in v.Vertices)
                                {
                                    if (vertplane.DistanceTo(p.Location) > 0)
                                    {
                                        voussoir.Add(p.Location);
                                    }
                                }
                                if (voussoir.Count > 0)
                                {
                                    selectedvoussoirs1.Add(v);
                                    side0Voussoirs.Add(v);
                                }
                            }
                        }

                        else
                        {
                            var vertplane = new Plane(springerLines[1].PointAtStart, springerLines[1].PointAtEnd, springerLines[1].PointAtStart + Plane.WorldXY.ZAxis * 1);
                            var prez = vertplane.ZAxis;
                            prez.Unitize();

                            var dLine = new Line(epicenter, line1.PointAt(0.5)).Direction;
                            dLine.Unitize();
                            if (prez * dLine < 0) vertplane.Flip();
                            List<Brep> flippedRow = new List<Brep>(orderedRowVoussoirs);
                            flippedRow.Reverse();

                            foreach (Brep v in flippedRow)
                            {
                                List<Point3d> voussoir = new List<Point3d>();
                                foreach (var p in v.Vertices)
                                {
                                    if (vertplane.DistanceTo(p.Location) > 0)
                                    {
                                        voussoir.Add(p.Location);
                                    }
                                }
                                if (voussoir.Count > 0)
                                {
                                    selectedvoussoirs2.Add(v);
                                    side1Voussoirs.Add(v);

                                }
                            }
                        }
                    }

                    side0.Add(selectedvoussoirs1.Count);
                    side1.Add(selectedvoussoirs2.Count);
                }
                var side0Max = side0.Max();
                var side1Max = side1.Max();                               

                double maxZ1 = double.MinValue;
                foreach (var p in side0Voussoirs)
                {
                    foreach (var t in p.Vertices)
                    {
                        if (t.Location.Z > maxZ1)
                        {
                            maxZ1 = t.Location.Z;
                        }
                    }
                }

                double maxZ2 = double.MinValue;
                foreach (var p in side1Voussoirs)
                {
                    foreach (var t in p.Vertices)
                    {
                        if (t.Location.Z > maxZ2)
                        {
                            maxZ2 = t.Location.Z;
                        }
                    }
                }

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

                    int q = side0Max; // remove first x
                    int w = side1Max; // remove last y

                    if (culledvoussoirs.Count >= q + w)
                    {
                        culledvoussoirs.RemoveRange(culledvoussoirs.Count - w, w); // remove last y
                        culledvoussoirs.RemoveRange(0, q);               // remove first x
                    }
                    outvoussoirs.AppendRange(culledvoussoirs.Select(v => new GH_Brep(v)), path);

                    List<Brep> selectedvoussoirs1 = new List<Brep>();
                    List<Brep> selectedvoussoirs2 = new List<Brep>();

                    for (int j = 0; j < 2; j++)
                    {
                        if (j == 0)
                        {                            
                            for (int k = 0; k < side0Max; k++)
                            {
                                var voussoir = orderedRowVoussoirs[k];                                
                                selectedvoussoirs1.Add(voussoir);                                
                            }
                        }

                        else
                        {                            
                            List<Brep> flippedRow = new List<Brep>(orderedRowVoussoirs);
                            flippedRow.Reverse();

                            for (int k = 0; k < side1Max; k++)
                            {
                                var voussoir = flippedRow[k];
                                selectedvoussoirs2.Add(voussoir);
                            }
                        }
                    }

                    b1.AppendRange(selectedvoussoirs1.Select(x => new GH_Brep(x)));
                    b3.AppendRange(selectedvoussoirs2.Select(x => new GH_Brep(x)));

                    //Intrados - Extrados
                    List<Brep> intradosside0 = new List<Brep>();
                    List<Brep> extradosside0 = new List<Brep>();
                    List<Brep> intradosside1 = new List<Brep>();
                    List<Brep> extradosside1 = new List<Brep>();
                    for (int l = 0; l < 2; l++)
                    {
                        if (l == 0)
                        {
                            foreach (Brep v in selectedvoussoirs1)
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
                                    .Select((f, p) => new { Face = f, Product = products[p] })
                                    .OrderBy(x => x.Product)
                                    .Select(x => x.Face)
                                    .ToList();

                                intradosside0.Add(faces[0]);
                                extradosside0.Add(faces[1]);
                            }
                        }
                        else
                        {
                            foreach (Brep v in selectedvoussoirs2)
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
                                    .Select((f, p) => new { Face = f, Product = products[p] })
                                    .OrderBy(x => x.Product)
                                    .Select(x => x.Face)
                                    .ToList();

                                intradosside1.Add(faces[0]);
                                extradosside1.Add(faces[1]);
                            }
                        }
                    }
                    //log.AppendRange(extrados.Select(x => new GH_Brep(x)), path);
                    b2.AppendRange(extradosside0.Select(x => new GH_Brep(x)));

                    for (int j = 0; j < springerLines.Count; j++)
                    {

                        if (j == 0)
                        {
                            List<Point3d> springerPolylinePoints1 = new List<Point3d>();
                            List<Point3d> springerPolylinePoints2 = new List<Point3d>();                                                         

                            for (int k = 0; k < 2; k++)
                            {
                                var intpointA = Intersection.CurvePlane(springerLines[1], transPlanes[i + k], RhinoMath.ZeroTolerance)[0].PointA;
                             
                                var transformExtrados = extradosside0.Last();
                                //b2.Append(new GH_Brep(transformExtrados));

                                List<Point3d> springerPLinePoints = new List<Point3d>();

                                //Pt 1
                                var pt1 = Intersection.CurvePlane(springerLines[0], transPlanes[i + k], RhinoMath.ZeroTolerance)[0].PointA;
                                //p1.Append(new GH_Point(pt1), path);

                                //Pt 2
                                //Pt 2
                                List<Point3d> ptsA = new List<Point3d>();

                                foreach (var v in intradosside0)
                                {
                                    var pt = v.Vertices
                                        .Select(vx => vx.Location)
                                        .OrderBy(p => Math.Abs(transPlanes[i + k].DistanceTo(p)))
                                        .Take(2)
                                        .OrderBy(p => p.DistanceTo(pt1))
                                        .ToList();

                                    ptsA.AddRange(pt);
                                }
                                ptsA.RemoveAt(0);
                                //outspringers.Append(new GH_Point(pt2), path);
                                //p1.AppendRange(ptsA.Select(x => new GH_Point(x)), path);

                                //Pt 3
                                var pt3 = transformExtrados.Vertices
                                    .Select(vx => vx.Location)
                                    .OrderBy(p => Math.Abs(transPlanes[i + k].DistanceTo(p)))
                                    .Take(2)
                                    .OrderBy(p => p.DistanceTo(pt1))
                                    .Last();
                                //p1.Append(new GH_Point(pt3), path);

                                //Pt 4                                
                                var pt4 = new Point3d(pt1.X, pt1.Y, maxZ1);
                                p1.Add(pt4);

                                //Pt 5
                                var pt5 = pt4 + intLines[i + k] * width;
                                //p1.Append(new GH_Point(pt5), path);

                                //Pt 6
                                var pt6 = pt1 + intLines[i + k] * width;
                                //p1.Append(new GH_Point(pt6), path);

                                if (k == 0)
                                {
                                    springerPolylinePoints1.Add(pt1);
                                    springerPolylinePoints1.AddRange(ptsA);
                                    springerPolylinePoints1.Add(pt3);
                                    springerPolylinePoints1.Add(pt4);
                                    springerPolylinePoints1.Add(pt5);
                                    springerPolylinePoints1.Add(pt6);
                                    springerPolylinePoints1.Add(pt1);
                                }
                                else
                                {
                                    springerPolylinePoints2.Add(pt1);
                                    springerPolylinePoints2.AddRange(ptsA);
                                    springerPolylinePoints2.Add(pt3);
                                    springerPolylinePoints2.Add(pt4);
                                    springerPolylinePoints2.Add(pt5);
                                    springerPolylinePoints2.Add(pt6);
                                    springerPolylinePoints2.Add(pt1);
                                }
                            }

                            // Assume springerPolylinePoints1 and springerPolylinePoints2 are your point lists
                           
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
                            var springerBody = Brep.JoinBreps(openSpringer, 0.01);

                            Plane sppPlane;
                            Plane.FitPlaneToPoints(springerPolylinePoints1, out sppPlane);
                            var springerProjectedPts = springerPolylinePoints1
                                .Select(p => sppPlane.ClosestPoint(p))
                                .ToList();

                            Plane sppPlane1;
                            Plane.FitPlaneToPoints(springerPolylinePoints2, out sppPlane1);
                            var springerProjectedPts1 = springerPolylinePoints2
                                .Select(p => sppPlane1.ClosestPoint(p))
                                .ToList();
                            var pol1 = new Polyline(springerProjectedPts).ToNurbsCurve();
                            c1.Append(new GH_Curve(pol1));
                            var pol2 = new Polyline(springerProjectedPts1).ToNurbsCurve();
                            c3.Append(new GH_Curve(pol2));

                            // Create planar caps for top and bottom
                            var face1 = Brep.CreatePlanarBreps(pol1, 0.001);
                            var face2 = Brep.CreatePlanarBreps(pol2, 0.001);

                            // Join the loft body with the planar caps
                            List<Brep> springerParts = new List<Brep>();
                            springerParts.Add(face1[0]);
                            springerParts.Add(springerBody[0]);
                            springerParts.Add(face2[0]);

                            Brep[] joinedSpringer = Brep.JoinBreps(springerParts, 0.01);

                            outspringers.Append(new GH_Brep(joinedSpringer[0]), path);
                        }

                        else
                        {
                            List<Point3d> springerPolylinePoints1 = new List<Point3d>();
                            List<Point3d> springerPolylinePoints2 = new List<Point3d>();

                            for (int k = 0; k < 2; k++)
                            {
                                var intpointA = Intersection.CurvePlane(springerLines[0], transPlanes[i + k], RhinoMath.ZeroTolerance)[0].PointA;

                                var transformExtrados = extradosside1.Last();
                                //b2.Append(new GH_Brep(transformExtrados));

                                List<Point3d> springerPLinePoints = new List<Point3d>();

                                //Pt 1
                                var pt1 = Intersection.CurvePlane(springerLines[1], transPlanes[i + k], RhinoMath.ZeroTolerance)[0].PointA;
                                //outspringers.Append(new GH_Point(pt1), path);

                                //Pt 2
                                List<Point3d> ptsA = new List<Point3d>();

                                foreach (var v in intradosside1)
                                {
                                    var pt = v.Vertices
                                        .Select(vx => vx.Location)
                                        .OrderBy(p => Math.Abs(transPlanes[i + k].DistanceTo(p)))
                                        .Take(2)
                                        .OrderBy(p => p.DistanceTo(pt1))
                                        .ToList();
                                                                        
                                    ptsA.AddRange(pt);
                                }
                                ptsA.RemoveAt(0);
                                //outspringers.Append(new GH_Point(pt2), path);
                                //p1.AppendRange(ptsA.Select(x => new GH_Point(x)), path);

                                //Pt 3
                                var pt3 = transformExtrados.Vertices
                                    .Select(vx => vx.Location)
                                    .OrderBy(p => Math.Abs(transPlanes[i + k].DistanceTo(p)))
                                    .Take(2)
                                    .OrderBy(p => p.DistanceTo(pt1))
                                    .Last();
                                //p1.Append(new GH_Point(pt3), path);

                                //Pt 4                                
                                var pt4 = new Point3d(pt1.X, pt1.Y, maxZ2);
                                p2.Add(pt4);

                                //Pt 5
                                var pt5 = pt4 + -intLines[i + k] * width;
                                
                                //Pt 6
                                var pt6 = pt1 + -intLines[i + k] * width;

                                if (k == 0)
                                {
                                    springerPolylinePoints1.Add(pt1);
                                    springerPolylinePoints1.AddRange(ptsA);
                                    springerPolylinePoints1.Add(pt3);
                                    springerPolylinePoints1.Add(pt4);
                                    springerPolylinePoints1.Add(pt5);
                                    springerPolylinePoints1.Add(pt6);
                                    springerPolylinePoints1.Add(pt1);
                                }
                                else
                                {
                                    springerPolylinePoints2.Add(pt1);
                                    springerPolylinePoints2.AddRange(ptsA);
                                    springerPolylinePoints2.Add(pt3);
                                    springerPolylinePoints2.Add(pt4);
                                    springerPolylinePoints2.Add(pt5);
                                    springerPolylinePoints2.Add(pt6);
                                    springerPolylinePoints2.Add(pt1);
                                }
                            }

                            // Assume springerPolylinePoints1 and springerPolylinePoints2 are your point lists

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

                            Plane sppPlane;
                            Plane.FitPlaneToPoints(springerPolylinePoints1, out sppPlane);
                            var springerProjectedPts = springerPolylinePoints1
                                .Select(p => sppPlane.ClosestPoint(p))
                                .ToList();

                            Plane sppPlane1;
                            Plane.FitPlaneToPoints(springerPolylinePoints2, out sppPlane1);
                            var springerProjectedPts1 = springerPolylinePoints2
                                .Select(p => sppPlane1.ClosestPoint(p))
                                .ToList();
                            var pol1 = new Polyline(springerProjectedPts).ToNurbsCurve();
                            c1.Append(new GH_Curve(pol1));
                            var pol2 = new Polyline(springerProjectedPts1).ToNurbsCurve();
                            c3.Append(new GH_Curve(pol2));

                            // Create planar caps for top and bottom
                            var face1 = Brep.CreatePlanarBreps(pol1, 0.001);
                            var face2 = Brep.CreatePlanarBreps(pol2, 0.001);

                            // Join the loft body with the planar caps
                            List<Brep> springerParts = new List<Brep>();
                            springerParts.Add(face1[0]);
                            springerParts.Add(springerBody[0]);
                            springerParts.Add(face2[0]);

                            Brep[] joinedSpringer = Brep.JoinBreps(springerParts, 0.01);

                            outspringers.Append(new GH_Brep(joinedSpringer[0]), path);
                        }
                    }
                }
                var p0l = new Line(p1.Last(), p1[0]);
                WL.Append(new GH_Curve(new LineCurve(p0l)), path);
                var p1l = new Line(p2.Last(), p2[0]);
                WL.Append(new GH_Curve(new LineCurve(p1l)), path);

            }
            DA.SetDataTree(0, outspringers);
            DA.SetDataTree(1, outvoussoirs);
            DA.SetDataTree(2, WL);
            //DA.SetDataTree(3, p1);
            //DA.SetDataTree(4, p2);
            //DA.SetDataTree(5, p3);
            //DA.SetDataTree(6, c1);
            //DA.SetDataTree(7, c2);
            //DA.SetDataTree(8, c3);
            //DA.SetDataTree(9, b1);
            //DA.SetDataTree(10, b2);
            //DA.SetDataTree(11, b3);

        }
    }
}
