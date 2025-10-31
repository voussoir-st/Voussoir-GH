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
        public static void SplitTreeAt<T>(
            GH_Structure<T> inputTree, 
            int splitIndex,
            out GH_Structure<T> treeA, 
            out GH_Structure<T> treeB)
            where T : IGH_Goo
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
            pManager.AddBrepParameter("Voussoirs", "V", "Non transformed voussoirs", GH_ParamAccess.item);
            //pManager.AddPointParameter("Log", "L", "All messages generated during execution", GH_ParamAccess.tree);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            //==============================
            // Inputs
            //==============================
            List<Curve> spLines = new List<Curve>();
            if (!DA.GetDataList(0, spLines)) return;
            GH_Structure<GH_Brep> remappedVoussoirs = new GH_Structure<GH_Brep>();
            if (!DA.GetDataTree(1, out remappedVoussoirs)) return;
            List<GH_Plane> tPlanes = new List<GH_Plane>();
            if (!DA.GetDataList(2, tPlanes)) return;
            double spWidth = 0.3;
            if (!DA.GetData(3, ref spWidth)) return;

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

            List<Point3d> avgPoints = new List<Point3d>();

            foreach (var branch in remappedVoussoirs.Branches)
            {
                foreach (var ghBrep in branch)
                {
                    Brep brep = ghBrep.Value;
                    if (brep == null || brep.Vertices.Count == 0)
                        continue;

                    double sumX = 0, sumY = 0, sumZ = 0;
                    foreach (var v in brep.Vertices)
                    {
                        sumX += v.Location.X;
                        sumY += v.Location.Y;
                        sumZ += v.Location.Z;
                    }

                    int n = brep.Vertices.Count;
                    avgPoints.Add(new Point3d(sumX / n, sumY / n, sumZ / n));
                }
            }

            List<GeometryBase> geom = avgPoints.Select(p => new Point(p)).Cast<GeometryBase>().ToList();
            double tolerance = RhinoMath.ZeroTolerance;           

            // Create the patch surface
            Brep patch = Brep.CreatePatch(
                geom,
                null,
                tolerance
            );
            if (patch != null && patch.Faces.Count > 0)
            {
                PolylineUtils.AlignNormalToWorldZ(patch.Faces[0].UnderlyingSurface());
            }

            //==============================
            // Initialize output structures
            //==============================
            var almostSpringers = new GH_Structure<GH_Brep>();
            var restVoussoirs = new GH_Structure<GH_Brep>();


            //==============================
            // Main loop over springer lines
            //==============================
            GH_Structure<GH_Brep> springerVoussoirs = new GH_Structure<GH_Brep>();
            GH_Structure<GH_Brep> extradosTree = new GH_Structure<GH_Brep>();
            GH_Structure<GH_Brep> intradosTree = new GH_Structure<GH_Brep>();

            //foreach (var line in spLines)
            for (int ispringer = 0; ispringer < spLines.Count; ispringer++)
            {
                var line = spLines[ispringer];
                
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

                var springerLinesLocal = new GH_Structure<GH_Curve>();

                for (int i = 0; i < springerLines.Count - 1; i++)
                {
                    var lines = new List<GH_Curve> { new GH_Curve(springerLines[i].ToNurbsCurve()), new GH_Curve(springerLines[i + 1].ToNurbsCurve()) };
                    springerLinesLocal.AppendRange(lines, new GH_Path(i));
                }

                // Extrude line to create vertical reference
                var extrusionBrep = Brep.CreateFromSurface(Surface.CreateExtrusion(line, Vector3d.ZAxis * 10));

                // Intersect with voussoirs
                var intersectedVoussoirs = new GH_Structure<GH_Brep>();
                foreach (var branchIdx in Enumerable.Range(0, remappedVoussoirs.Branches.Count))
                {
                    var branch = remappedVoussoirs.Branches[branchIdx];
                    foreach (var ghBrep in branch)
                    {
                        var voussoir = ghBrep.Value;
                        var points = voussoir.Vertices;
                        var truepoints = new List<Point3d>();
                        foreach (var point in points)
                        {
                            if (spPlane.DistanceTo(point.Location) > 0) truepoints.Add(point.Location);
                        }
                        if (truepoints.Count > 0)
                        {
                            intersectedVoussoirs.Append(ghBrep, remappedVoussoirs.Paths[branchIdx]);                           
                        }
                    }
                }

                int maxCount = intersectedVoussoirs.Branches.Max(b => b.Count);
                GH_Structure<GH_Brep> spVoussoirs = new GH_Structure<GH_Brep>();

                for (int i = 0; i < remappedVoussoirs.Branches.Count; i++)
                {
                    var rowVoussoirs = remappedVoussoirs.Branches[i];

                    if (ispringer == 0)
                    {
                        for (int j = 0; j < maxCount; j++)
                        {
                            var brep = rowVoussoirs[j];
                            var sp = brep.Value;
                            springerVoussoirs.Append(new GH_Brep(sp), new GH_Path(i));
                            spVoussoirs.Append(new GH_Brep(sp), new GH_Path(i));
                        }
                    }
                    if (ispringer == 1)
                    {
                        for (int j = 0; j < maxCount; j++)
                        {
                            var brep = rowVoussoirs[rowVoussoirs.Count - 1 - j];
                            var sp = brep.Value;
                            springerVoussoirs.Append(new GH_Brep(sp), new GH_Path(i));
                            spVoussoirs.Append(new GH_Brep(sp), new GH_Path(i));
                        }
                    }
                }



                for (int i = 0; i < spVoussoirs.PathCount; i++)
                {
                    GH_Path path = spVoussoirs.get_Path(i);
                    var voussoirBranch = spVoussoirs.Branches[i];

                    foreach (var vouss in voussoirBranch)
                    {
                        Brep voussoir = vouss.Value;

                        // 1️⃣ Compute centroid (average of vertex positions)
                        var verts = voussoir.Vertices;

                        Point3d center = Point3d.Origin;
                        foreach (var c in verts)
                            center += c.Location;
                        center /= verts.Count;

                        // 2️⃣ Find closest point on patch
                        double u, v;
                        patch.Faces[0].ClosestPoint(center, out u, out v);

                        Point3d pt;
                        Vector3d[] derivs;
                        var normal = patch.Faces[0].NormalAt(u, v);
                        normal.Unitize();

                        // 3️⃣ Compare each face normal with patch normal
                        double maxDot = double.MinValue;
                        double minDot = double.MaxValue;
                        BrepFace extradosFace = null;
                        BrepFace intradosFace = null;

                        foreach (var face in voussoir.Faces)
                        {
                            double a, b;
                            face.ClosestPoint(center, out a, out b);
                            var faceNormal = face.NormalAt(a, b);
                            faceNormal.Unitize();
                            double product = normal * faceNormal;

                            if (product > maxDot)
                            {
                                maxDot = product;
                                extradosFace = face;
                            }
                            if (product < minDot)
                            {
                                minDot = product;
                                intradosFace = face;
                            }
                        }

                        // 4️⃣ Append to trees with matching structure
                        if (extradosFace != null)
                            extradosTree.Append(new GH_Brep(extradosFace.DuplicateFace(false)), path);

                        if (intradosFace != null)
                            intradosTree.Append(new GH_Brep(intradosFace.DuplicateFace(false)), path);
                    }
                }


                // Compute Z max
                double maxZ = double.MinValue;

                foreach (var branch in spVoussoirs.Branches)
                {
                    foreach (var ghBrep in branch)
                    {
                        var brep = ghBrep.Value;
                        foreach (var v in brep.Vertices)
                        {
                            double z = v.Location.Z;
                            if (z > maxZ)
                                maxZ = z;
                        }
                    }
                }

                for (int i = 0; i < remappedVoussoirs.Branches.Count; i++)
                {
                    for (int j = 0; j < 2; j++)
                    {
                        List<Point3d> spPoints = new List<Point3d>();

                        //Bottom out
                        Point3d pt1 = Point3d.Unset; 
                        pt1 = springerLines[i + j].PointAt(1);
                        spPoints.Add(pt1);

                        //Bottom in
                        Point3d pt2 = Point3d.Unset; 
                        pt2 = springerLines[i + j].PointAt(0);
                        spPoints.Add(pt2);

                        //Intrados * vertical plane
                        Point3d pt3 = Point3d.Unset;

                        //or
                        //Intrados * XYPlane
                        Point3d pt4 = Point3d.Unset;

                        //Intrados middle points
                        List<Point3d> ptsA = new List<Point3d>();

                        //Voussoir highest
                        Point3d pt5 = Point3d.Unset;

                        //Top in
                        Point3d pt6 = Point3d.Unset; 
                        pt6 = springerLines[i + j].PointAt(0) + Vector3d.ZAxis * maxZ;
                        spPoints.Add(pt6);

                        //Top out
                        Point3d pt7 = Point3d.Unset; 
                        pt7 = springerLines[i + j].PointAt(1) + Vector3d.ZAxis * maxZ;
                        spPoints.Add(pt7);

                    }
                }
            }
            GH_Structure<GH_Brep> trimspringerVoussoirs = TreeUtils.TrimTreeDepth(springerVoussoirs);
            
            //==============================
            // Output
            //==============================
            DA.SetDataTree(0, intradosTree);
            //DA.SetData(1, patch);
            //DA.SetDataTree(2, voussoirPoints);
        }
    }  
}
