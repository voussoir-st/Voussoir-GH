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
            //==============================
            // Inputs
            //==============================

            // Get input 0: list of springer base lines (probably defining arch boundaries or divisions)
            Brep surface = null;
            if (!DA.GetData(0, ref surface)) return; // If missing, stop execution

            // Get input 0: list of springer base lines (probably defining arch boundaries or divisions)
            List<Curve> spLines = new List<Curve>();
            if (!DA.GetDataList(1, spLines)) return; // If missing, stop execution

            // Get input 1: data tree containing remapped voussoirs (stone blocks as Breps)
            GH_Structure<GH_Brep> remappedVoussoirs = new GH_Structure<GH_Brep>();
            if (!DA.GetDataTree(2, out remappedVoussoirs)) return;

            // Get input 2: list of target planes, probably defining horizontal layers or construction levels
            List<GH_Plane> tPlanes = new List<GH_Plane>();
            if (!DA.GetDataList(3, tPlanes)) return;

            // Get input 3: width of the springer section (initial default = 0.3 units)
            double spWidth = 0.3;
            if (!DA.GetData(4, ref spWidth)) return;

            //==============================
            // Extract extrados faces from voussoirs
            //==============================
            // This section loops through all voussoirs (Brep stones) and extracts what is assumed
            // to be the “extrados” face (outer upper face) from each, creating a new Brep collection.

            var extrados = new GH_Structure<GH_Brep>();
            foreach (var b in remappedVoussoirs.AllData(true).OfType<GH_Brep>())
            {
                // Ensure the Brep has multiple faces before selecting
                if (b.Value.Faces.Count > 2)
                {
                    // Duplicate the face at index 1 (assumed to be extrados)
                    var extrBrep = b.Value.Faces[1].DuplicateFace(false);

                    // Store it in a new Grasshopper data structure (path 0 since no hierarchy is needed yet)
                    extrados.Append(new GH_Brep(extrBrep), new GH_Path(0));
                }
            }

            //==============================
            // Compute closest points & planar distance
            //==============================
            // This computes the minimal 2D distance between each springer line and the extrados surfaces,
            // then uses the maximum found distance to adjust the “spWidth” if it’s too small.

            double planarDistance = 0;
            foreach (var sLine in spLines)
            {
                // Create a flattened version of the line start point (Z=0) for 2D comparison
                Point3d startXY = new Point3d(sLine.PointAtStart.X, sLine.PointAtStart.Y, 0);
                double minDist = double.MaxValue;

                // For each extrados surface, compute the closest point to the start of this line
                foreach (var ghBrep in extrados.AllData(true).OfType<GH_Brep>())
                {
                    // Find the closest point on the Brep to the line’s start
                    if (ghBrep.Value.ClosestPoint(sLine.PointAtStart, out Point3d pt, out _, out _, out _, double.MaxValue, out _))
                    {
                        // Compare only the XY-plane distance
                        double distXY = new Point3d(pt.X, pt.Y, 0).DistanceTo(startXY);
                        minDist = Math.Min(minDist, distXY); // Keep smallest found
                    }
                }

                // Keep the largest minimal distance among all lines (the worst case)
                planarDistance = Math.Max(planarDistance, minDist);
            }

            // If the initial width is smaller than the found planar distance, expand it
            if (spWidth < planarDistance) spWidth = planarDistance;

            //==============================
            // Compute centroid of all vertices
            //==============================
            // This section computes the 3D centroid of all voussoir vertices in the entire structure,
            // to use as a global spatial reference for later operations (e.g., plane orientation).

            var allVerts = remappedVoussoirs.AllData(true).OfType<GH_Brep>()
                .Where(b => b.Value != null)
                .SelectMany(b => b.Value.Vertices.Select(v => v.Location))
                .ToList();

            // Compute average X, Y, Z of all vertices
            Point3d centroid = new Point3d(
                allVerts.Average(p => p.X),
                allVerts.Average(p => p.Y),
                allVerts.Average(p => p.Z)
            );

            //==============================
            // Compute average point for each voussoir
            //==============================
            // For each individual voussoir Brep, compute its own vertex average point
            // (used to build a cloud of representative points later).

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
                    avgPoints.Add(new Point3d(sumX / n, sumY / n, sumZ / n)); // Local average per voussoir
                }
            }

            // Create a Brep patch that fits the input points
            Brep patch = surface;

            if (patch != null && patch.Faces.Count > 0)
            {
                var face = patch.Faces[0];
                var surf = face.UnderlyingSurface();
                var du = surf.Domain(0);
                var dv = surf.Domain(1);
                double u = (du.T0 + du.T1) / 2.0;
                double v = (dv.T0 + dv.T1) / 2.0;
                var normal = surf.NormalAt(u, v);

                if (normal.Z < 0)
                {
                    // Flip entire Brep orientation
                    patch.Flip(); // ✅ this method exists in Brep
                }
            }

            //==============================
            // Initialize output structures
            //==============================
            // These GH_Structures will store categorized voussoirs (springers and the rest)
            var springers = new GH_Structure<GH_Brep>();
            var restVoussoirs = new GH_Structure<GH_Brep>();

            //==============================
            // Main loop over springer lines
            //==============================
            
            GH_Structure<GH_Brep> springerVoussoirs = new GH_Structure<GH_Brep>();
            var ptpt = new GH_Structure<GH_Curve>(); 
            var lcounter = 0; // just a debug counter to track iterations

            var removed0 = new GH_Structure<GH_Brep>();
            var removed1 = new GH_Structure<GH_Brep>();

            var wallLinesPoints1 = new List<Point3d>();
            var wallLinesPoints2 = new List<Point3d>();

            // Loop through each springer line
            for (int ispringer = 0; ispringer < spLines.Count; ispringer++)
            {
                //Debug.WriteLine($"\n\nlcounter: " + lcounter);  // Debug info for console tracking
                var line = spLines[ispringer];                   // Current springer line

                //==============================
                // Compute base surfaces (springerSurfaces)
                //==============================
                var spPlane = new Plane(line.PointAtStart, line.PointAtEnd - line.PointAtStart, Vector3d.ZAxis);

                // Flip the plane if it faces away from the overall geometry centroid
                if (spPlane.DistanceTo(centroid) > 0) spPlane.Flip();

                // Original line (used as base)
                Line sLine = new Line(line.PointAtStart, line.PointAtEnd);

                // Offset line along Z-axis (thickness or projection upward)
                var offsetLine = new Line(line.PointAtStart + spPlane.ZAxis * spWidth, line.PointAtEnd + spPlane.ZAxis * spWidth);
                offsetLine.Extend(20, 20); // Extend both ends for safety (ensures intersection with planes)

                var basePlane = new Plane(line.PointAtStart, line.PointAtEnd, offsetLine.PointAt(0));
                var baseZ = basePlane.ZAxis;
                baseZ.Unitize();

                if (baseZ.Z < 0)
                {
                    basePlane.Flip();
                }

                // Container for computed springer lines between tPlanes
                var springerLines = new List<Line>();

                //==============================
                // Intersect springer lines with target planes
                //==============================
                foreach (var p in tPlanes)
                {
                    double tA, tB;
                    bool hasA = Intersection.LinePlane(sLine, p.Value, out tA);          // intersection with base line
                    bool hasB = Intersection.LinePlane(offsetLine, p.Value, out tB);     // intersection with offset line

                    // If both intersections are valid, create a line segment between the two points
                    if (hasA && hasB)
                    {
                        springerLines.Add(new Line(sLine.PointAt(tA), offsetLine.PointAt(tB)));
                    }
                }

                //==============================
                // Create a hierarchical structure of local springer lines
                //==============================
                var springerLinesLocal = new GH_Structure<GH_Curve>();

                for (int i = 0; i < springerLines.Count - 1; i++)
                {
                    // Create curve pairs between consecutive lines for further geometry generation
                    var lines = new List<GH_Curve>
                    {
                        new GH_Curve(springerLines[i].ToNurbsCurve()),
                        new GH_Curve(springerLines[i + 1].ToNurbsCurve())
                    };
                    springerLinesLocal.AppendRange(lines, new GH_Path(i));  // Store with a unique path index
                }

                //==============================
                // Extrude line to create vertical reference surface
                //==============================
                var extrusionBrep = Brep.CreateFromSurface(Surface.CreateExtrusion(line, Vector3d.ZAxis * 10));

                //==============================
                // Intersect this extrusion with all voussoirs
                //==============================
                var intersectedVoussoirs = new GH_Structure<GH_Brep>();
                foreach (var branchIdx in Enumerable.Range(0, remappedVoussoirs.Branches.Count))
                {
                    var branch = remappedVoussoirs.Branches[branchIdx];
                    foreach (var ghBrep in branch)
                    {
                        var voussoir = ghBrep.Value;
                        var points = voussoir.Vertices;
                        var truepoints = new List<Point3d>();

                        // Keep vertices that are above the current springer plane
                        foreach (var point in points)
                        {
                            if (spPlane.DistanceTo(point.Location) >= 0 || spPlane.DistanceTo(point.Location) >= -0.01)
                                truepoints.Add(point.Location);
                        }

                        // If the voussoir intersects (some vertices above the plane), add it to structure
                        if (truepoints.Count > 0)
                        {
                            intersectedVoussoirs.Append(ghBrep, remappedVoussoirs.Paths[branchIdx]);
                        }
                    }
                }

                //==============================
                // Determine maximum number of voussoirs per branch
                //==============================
                int maxCount = intersectedVoussoirs.Branches.Max(b => b.Count);

                // Temporary structure to hold selected voussoirs for this springer
                GH_Structure<GH_Brep> spVoussoirs = new GH_Structure<GH_Brep>();

                //==============================
                // Handle orientation of voussoirs per springer line
                //==============================
                for (int i = 0; i < remappedVoussoirs.Branches.Count; i++)
                {
                    var rowVoussoirs = remappedVoussoirs.Branches[i];

                    // For the first springer, add voussoirs in normal order
                    if (ispringer == 0)
                    {
                        for (int j = 0; j < maxCount; j++)
                        {
                            var brep = rowVoussoirs[j];
                            var sp = brep.Value;
                            springerVoussoirs.Append(new GH_Brep(sp), new GH_Path(i));
                            spVoussoirs.Append(new GH_Brep(sp), new GH_Path(i));
                            removed0.Append(new GH_Brep(sp), new GH_Path(i));                            
                        }
                    }

                    // For the second springer, reverse order (mirroring logic)
                    if (ispringer == 1)
                    {
                        for (int j = 0; j < maxCount; j++)
                        {
                            var brep = rowVoussoirs[rowVoussoirs.Count - 1 - j];
                            var sp = brep.Value;
                            springerVoussoirs.Append(new GH_Brep(sp), new GH_Path(i));
                            spVoussoirs.Append(new GH_Brep(sp), new GH_Path(i));
                            removed1.Append(new GH_Brep(sp), new GH_Path(i));
                        }
                    }
                }

                //==============================
                // Create trees for extrados and intrados faces
                //==============================
                GH_Structure<GH_Brep> extradosTree = new GH_Structure<GH_Brep>();
                GH_Structure<GH_Brep> intradosTree = new GH_Structure<GH_Brep>();

                // Loop through all voussoirs to separate top and bottom surfaces
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

                        // 2️⃣ Find closest point on patch surface (global reference)
                        double u, v;
                        patch.Faces[0].ClosestPoint(center, out u, out v);

                        var normal = patch.Faces[0].NormalAt(u, v);
                        normal.Unitize();

                        // 3️⃣ Compare each face normal with patch normal
                        // The face most aligned (dot product max) → extrados
                        // The face most opposite (dot product min) → intrados
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

                        // 4️⃣ Append to trees (keeping same GH_Path for consistency)
                        if (extradosFace != null)
                            extradosTree.Append(new GH_Brep(extradosFace.DuplicateFace(false)), path);

                        if (intradosFace != null)
                            intradosTree.Append(new GH_Brep(intradosFace.DuplicateFace(false)), path);
                    }
                }

                //==============================
                // Compute maximum Z elevation among all voussoir vertices
                //==============================
                double maxZ = double.MinValue;

                foreach (var branch in spVoussoirs.Branches)
                {
                    foreach (var ghBrep in branch)
                    {
                        var brep = ghBrep.Value;
                        foreach (var v in brep.Vertices)
                        {
                            double z = v.Location.Z;
                            //Debug.WriteLine($"z: " + z);
                            if (z > maxZ)
                                maxZ = z;
                        }
                    }
                }
                //Debug.WriteLine($"\nmaxZ: " + maxZ);


                //==============================
                // Loop through each voussoir branch (per springer)
                //==============================
                for (int i = 0; i < remappedVoussoirs.Branches.Count; i++)
                {
                    // --- Collect intrados face vertices ---
                    var branchintrados = intradosTree.Branches[i];
                    var intradospoints1 = new List<Point3d>();
                    var intradospoints2 = new List<Point3d>();

                    foreach (var ghBrep in branchintrados)
                    {
                        var brep = ghBrep.Value;
                        var pointsA = new List<Point3d>();
                        var pointsB = new List<Point3d>();
                        //Debug.WriteLine($"brep vert: " + brep.Vertices.Count);

                        // Split vertices into two groups depending on which side of the tPlane they are
                        foreach (var v in brep.Vertices)
                        {
                            if (Math.Abs(tPlanes[i].Value.DistanceTo(v.Location)) < 1e-3)
                            {
                                // Vertex lies directly on the plane
                                pointsA.Add(v.Location);
                                //Debug.WriteLine($"z: " + Math.Abs(tPlanes[i].Value.DistanceTo(v.Location)));
                            }
                            else
                                pointsB.Add(v.Location);
                        }

                        // Sort both lists by Z to keep geometry consistent vertically
                        pointsA = pointsA.OrderBy(pt => pt.Z).ToList();
                        pointsB = pointsB.OrderBy(pt => pt.Z).ToList();

                        // Accumulate for later polyline construction
                        intradospoints1.AddRange(pointsA);
                        intradospoints2.AddRange(pointsB);
                    }

                    // --- Collect extrados face vertices ---
                    var branchextrados = extradosTree.Branches[i];
                    var extradospoints1 = new List<Point3d>();
                    var extradospoints2 = new List<Point3d>();

                    foreach (var ghBrep in branchextrados)
                    {
                        var brep = ghBrep.Value;
                        var pointsA = new List<Point3d>();
                        var pointsB = new List<Point3d>();

                        foreach (var v in brep.Vertices)
                        {
                            if (Math.Abs(tPlanes[i].Value.DistanceTo(v.Location)) < 0.001)
                                pointsA.Add(v.Location);
                            else
                                pointsB.Add(v.Location);
                        }
                        //Debug.WriteLine($"pointsA: " + pointsA.Count);
                        extradospoints1.AddRange(pointsA);
                        extradospoints2.AddRange(pointsB);
                    }

                    //==============================
                    // Generate section polylines per springer
                    //==============================
                    var linecounter = 0;
                    var polylines = new List<Curve>();
                    var springerDefLines0 = new List<Line>();
                    var springerDefLines1 = new List<Line>();

                    for (int j = 0; j < 2; j++) // two sides per voussoir row
                    {
                        //Debug.WriteLine($"\nlinecounter: " + linecounter);
                        List<Point3d> spPoints = new List<Point3d>();

                        // --- Bottom outer point (lower edge of springer) ---
                        Point3d pt1 = springerLines[i + j].PointAt(1);
                        spPoints.Add(pt1);

                        // --- Bottom inner point ---
                        Point3d pt2 = springerLines[i + j].PointAt(0);
                        spPoints.Add(pt2);

                        //// --- Intrados intersection with vertical extrusion (currently only computed, not used) ---
                        //Point3d ptB = Point3d.Unset;

                        //for (int a = 0; a < remappedVoussoirs.Branches.Count; a++)
                        //{
                        //    var bintrados = intradosTree.Branches[i];
                        //    foreach (var intrados in bintrados)
                        //    {
                        //        Surface srfA = intrados.Value.Surfaces[0];
                        //        Surface srfB = extrusionBrep.Surfaces[0];

                        //        Curve[] intersectionCurves;
                        //        Point3d[] intersectionPoints;

                        //        double tol = RhinoDoc.ActiveDoc.ModelAbsoluteTolerance;

                        //        bool success = Intersection.SurfaceSurface(srfA, srfB, tol, out intersectionCurves, out intersectionPoints);
                        //        // intersectionCurves and intersectionPoints could later define pt3
                        //    }
                        //}

                        //// Optional: Intrados * XYPlane intersection placeholder

                        // --- Intrados midpoints (used to trace interior profile) ---
                        List<Point3d> ptsA = new List<Point3d>();
                        ptsA = (j == 0) ? intradospoints1.Skip(1).ToList() : intradospoints2.Skip(1).ToList();


                        // --- Determine the highest point from extrados vertices ---
                        Point3d pt5 = Point3d.Unset;
                        if (j == 0)
                        {
                            var highest = Point3d.Unset;
                            var maxh = double.MinValue;
                            foreach (var p in extradospoints1)
                            {
                                if (p.Z > maxh)
                                {
                                    maxh = p.Z;
                                    highest = p;
                                }
                            }
                            pt5 = highest;
                        }
                        else
                        {
                            var highest = Point3d.Unset;
                            var maxh = double.MinValue;
                            foreach (var p in extradospoints2)
                            {
                                if (p.Z > maxh)
                                {
                                    maxh = p.Z;
                                    highest = p;
                                }
                            }
                            pt5 = highest;
                        }

                        // --- Top inner & outer points (raised by maxZ) ---
                        var pt6x = springerLines[i + j].PointAt(0).X;
                        var pt6y = springerLines[i + j].PointAt(0).Y;
                        Point3d pt6 = new Point3d(pt6x, pt6y, maxZ);
                        spPoints.Add(pt6);
                        if (ispringer == 0)
                        {
                            wallLinesPoints1.Add(pt6);
                        }
                        else
                        {
                            wallLinesPoints2.Add(pt6);
                        }

                        var pt7x = springerLines[i + j].PointAt(1).X;
                        var pt7y = springerLines[i + j].PointAt(1).Y;
                        Point3d pt7 = new Point3d(pt7x, pt7y, maxZ);
                        spPoints.Add(pt7);

                        //==============================
                        // Build polyline through ordered points
                        //==============================
                        var sPoints = new List<GH_Point>();
                        sPoints.Add(new GH_Point(pt1));
                        sPoints.Add(new GH_Point(pt2));
                        sPoints.AddRange(ptsA.Select(p => new GH_Point(p)));
                        sPoints.Add(new GH_Point(pt5));
                        sPoints.Add(new GH_Point(pt6));
                        sPoints.Add(new GH_Point(pt7));
                        sPoints.Add(new GH_Point(pt1)); // close polyline loop

                        //Convert GH_Points to Rhino Points and then to PolylineCurve
                        var pts = sPoints.Select(p => p.Value).ToList();
                        var curve = new PolylineCurve(new Polyline(pts));
                        polylines.Add(curve);
                        //Debug.WriteLine($"polylines: " + polylines.Count);
                        linecounter++;
                    }

                    //==============================
                    // Store section polylines and create 3D geometry
                    //==============================
                    foreach (var poly in polylines)
                    {
                        ptpt.Append(new GH_Curve(poly), new GH_Path(i));
                        //Debug.WriteLine($"ptpt: " + ptpt.PathCount);
                    }
                    //Debug.WriteLine($"polylines: " + polylines.Count);

                    // Loft between the two polylines to form the springer’s body
                    var sSpringer = TreeUtils.LoftBySegments(polylines[0], polylines[1]);

                    // Cap top and bottom faces with planar surfaces
                    var face1 = Brep.CreatePlanarBreps(polylines[0], RhinoMath.ZeroTolerance);
                    var face2 = Brep.CreatePlanarBreps(polylines[1], RhinoMath.ZeroTolerance);

                    // Join loft + faces into a closed brep
                    var springerlist = new List<Brep>();
                    springerlist.Add(face1[0]);
                    springerlist.Add(sSpringer);
                    springerlist.Add(face2[0]);

                    var joinedSpringer = Brep.JoinBreps(springerlist, RhinoMath.ZeroTolerance);

                    // Append to final springers output
                    springers.Append(new GH_Brep(joinedSpringer[0]), new GH_Path(i));
                }
                lcounter++;
            }
            var pta = wallLinesPoints1.First();
            var ptb = wallLinesPoints1.Last();
            var ptc = wallLinesPoints2.First();
            var ptd = wallLinesPoints2.Last();

            var topLines = new List<Line>();

            Line l1 = new Line(pta, ptb);
            if (l1.From.DistanceTo(spLines[0].PointAtStart) > l1.From.DistanceTo(spLines[0].PointAtEnd))
                l1.Flip();
            topLines.Add(l1);

            Line l2 = new Line(ptc, ptd);
            if (l2.From.DistanceTo(spLines[1].PointAtStart) > l2.From.DistanceTo(spLines[1].PointAtEnd))
                l2.Flip();
            topLines.Add(l2);


            //==============================
            // Separate the remaining voussoirs (not part of springers)
            //==============================
            for (int i = 0; i < remappedVoussoirs.Branches.Count; i++)
            {
                var path = remappedVoussoirs.Paths[i];
                var branch = remappedVoussoirs.Branches[i];

                // Determine how many to remove on each side
                int removeLeft = (i < removed0.Branches.Count) ? removed0.Branches[i].Count : 0;
                int removeRight = (i < removed1.Branches.Count) ? removed1.Branches[i].Count : 0;

                // Skip small branches (not enough voussoirs)
                if (branch.Count <= removeLeft + removeRight)
                    continue;

                // Keep only the middle voussoirs
                var kept = branch
                    .Skip(removeLeft)
                    .Take(branch.Count - removeLeft - removeRight)
                    .ToList();

                foreach (var brep in kept)
                    restVoussoirs.Append(brep, path);
            }

            //==============================
            // Output results to Grasshopper
            //==============================
            DA.SetDataTree(0, springers);      // 0 → complete springer geometry
            DA.SetDataTree(1, restVoussoirs);  // 1 → remaining voussoirs
            DA.SetDataList(2, topLines);
        }
    }  
}
