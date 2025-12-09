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
            //pManager.AddPointParameter("Log", "L", "All messages generated during execution", GH_ParamAccess.tree);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // --------------------------------------------
            // 1. READ INPUTS AS TREES/LISTS
            // --------------------------------------------

            // Vault surfaces as list
            GH_Structure<GH_Brep> vaultSurfaces = new GH_Structure<GH_Brep>();
            if (!DA.GetDataTree(0, out vaultSurfaces)) return;

            // Springer lines: tree
            GH_Structure<GH_Curve> spLinesTree;
            if (!DA.GetDataTree(1, out spLinesTree)) return;

            // Voussoirs: tree
            GH_Structure<GH_Brep> voussoirsFullTree;
            if (!DA.GetDataTree(2, out voussoirsFullTree)) return;

            // Transversal planes: tree
            GH_Structure<GH_Plane> planesTree;
            if (!DA.GetDataTree(3, out planesTree)) return;         
                
            // Outputs (multi-vault)
            GH_Structure<GH_Brep> outSpringers = new GH_Structure<GH_Brep>();
            GH_Structure<GH_Brep> outRest = new GH_Structure<GH_Brep>();


            // --------------------------------------------
            // 2. MULTI-VAULT LOOP
            // --------------------------------------------

            for (int v = 0; v < vaultSurfaces.Branches.Count; v++)
            {
                // --------------------------------------------
                // 2.1 Extract ONE vault's inputs
                // --------------------------------------------
           
                Brep thisVaultSurface = vaultSurfaces.Branches[v][0].Value;

                // Springer lines
                List<Curve> thisVaultSpringers =
                    spLinesTree.Branches[v].Select(c => c.Value).ToList();

                // Planes
                List<Plane> thisVaultPlanes =
                    planesTree.Branches[v].Select(p => p.Value).ToList();

                // Create structure for this vault - preserve original paths
                GH_Structure<GH_Brep> vousForVault = new GH_Structure<GH_Brep>();

                foreach (GH_Path p in voussoirsFullTree.Paths)
                {
                    // Check if this is a multi-vault structure by looking at path depth
                    bool isMultiVault = p.Indices.Length >= 3;

                    if (isMultiVault)
                    {
                        // Multi-vault case: {vault;A;B} - filter by vault index
                        if (p.Indices[0] != v) continue;
                    }
                    else
                    {
                        // Single vault case: {A;B} - all paths belong to vault 0
                        // Only process when v == 0, skip other iterations
                        if (v != 0) break; // Exit the vault loop entirely after first iteration
                    }

                    // Get the corresponding branch
                    var branchIList = voussoirsFullTree.get_Branch(p);
                    if (branchIList == null) continue;

                    // Convert IList → IEnumerable<GH_Brep>
                    var branch = branchIList.Cast<GH_Brep>();

                    // Add the branch using the ORIGINAL path
                    vousForVault.AppendRange(branch, p);
                }

                //Debug.WriteLine($"Vault {v}: vousForVault branches: " + vousForVault.Branches.Count);
                //if (vousForVault.Branches.Count > 0)
                //    Debug.WriteLine($"First branch count: " + vousForVault.Branches[0].Count);

                // --- Remap voussoirs: extract last two indices from original path ---
                var remappedVoussoirs = new GH_Structure<GH_Brep>();

                for (int branchIdx = 0; branchIdx < vousForVault.Paths.Count; branchIdx++)
                {
                    var oldPath = vousForVault.Paths[branchIdx];
                    var branch = vousForVault.Branches[branchIdx];

                    int A, B;
                    int pathLength = oldPath.Indices.Length;

                    // Extract the last two indices from the original path
                    if (pathLength >= 2)
                    {
                        A = oldPath.Indices[pathLength - 2];
                        B = oldPath.Indices[pathLength - 1];
                    }
                    else if (pathLength == 1)
                    {
                        A = 0;
                        B = oldPath.Indices[0];
                    }
                    else
                    {
                        A = 0;
                        B = 0;
                    }

                    //Debug.WriteLine($"Original path: {oldPath} -> Remapped to: {{B={B}}}(A={A})");

                    // Insert with remapped structure: {B}(A)
                    for (int i = 0; i < branch.Count; i++)
                    {
                        remappedVoussoirs.Insert(new GH_Brep(branch[i]), new GH_Path(B), A);
                    }
                }

                Debug.WriteLine($"Remapped - Paths Count: " + remappedVoussoirs.Branches.Count);
                Debug.WriteLine($"Remapped - First branch Count: " + (remappedVoussoirs.Branches.Count > 0 ? remappedVoussoirs.Branches[0].Count : 0));
                
                //foreach (GH_Path p in remappedVoussoirs.Paths)
                //{
                //     Debug.WriteLine($"remappedVoussoirs Path: " + p + ", p: " + remappedVoussoirs.get_Branch(p).Count);
                //}

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
                foreach (var sLine in thisVaultSpringers)
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
                var spWidth = 0.3;
                if (spWidth < planarDistance) spWidth = planarDistance;

                //==============================
                // Compute centroid of all vertices
                //==============================
                // This section computes the 3D centroid of all voussoir vertices in the entire structure,
                // to use as a global spatial reference for later operations (e.g., plane orientation).

                var allVerts = remappedVoussoirs.AllData(true).OfType<GH_Brep>()
                    .Where(b => b.Value != null)
                    .SelectMany(b => b.Value.Vertices.Select(z => z.Location))
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
                        foreach (var l in brep.Vertices)
                        {
                            sumX += l.Location.X;
                            sumY += l.Location.Y;
                            sumZ += l.Location.Z;
                        }

                        int n = brep.Vertices.Count;
                        avgPoints.Add(new Point3d(sumX / n, sumY / n, sumZ / n)); // Local average per voussoir
                    }
                }

                // Create a Brep patch that fits the input points
                Brep patch = thisVaultSurface;

                if (patch != null && patch.Faces.Count > 0)
                {
                    var face = patch.Faces[0];
                    var surf = face.UnderlyingSurface();
                    var du = surf.Domain(0);
                    var dv = surf.Domain(1);
                    double i = (du.T0 + du.T1) / 2.0;
                    double j = (dv.T0 + dv.T1) / 2.0;
                    var normal = surf.NormalAt(i, j);

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

                // Loop through each springer line
                for (int ispringer = 0; ispringer < thisVaultSpringers.Count; ispringer++)
                {
                    Debug.WriteLine($"\n\nispringer: " + ispringer);  // Debug info for console tracking
                    var line = thisVaultSpringers[ispringer];                   // Current springer line

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
                    //Debug.WriteLine($"tPlanes: " + thisVaultPlanes.Count);
                    foreach (var p in thisVaultPlanes)
                    {
                        
                        double tA, tB;
                        bool hasA = Intersection.LinePlane(sLine, p, out tA);          // intersection with base line
                        bool hasB = Intersection.LinePlane(offsetLine, p, out tB);     // intersection with offset line

                        // If both intersections are valid, create a line segment between the two points
                        if (hasA && hasB)
                        {
                            springerLines.Add(new Line(sLine.PointAt(tA), offsetLine.PointAt(tB)));
                        }
                    }
                    //Debug.WriteLine($"Springerlines: " + springerLines.Count);

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
                    // Determine maximum number of voussoirs per branch
                    //==============================
                    int maxCount = 1;

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

                    foreach (GH_Path p in spVoussoirs.Paths)
                        Debug.WriteLine($"spVoussoirs: " + p);

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
                            double n, m;
                            patch.Faces[0].ClosestPoint(center, out n, out m);

                            var normal = patch.Faces[0].NormalAt(n, m);
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
                            foreach (var n in brep.Vertices)
                            {
                                double z = n.Location.Z;
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
                    var auxLinePoints = new List<Point3d>();
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
                            foreach (var n in brep.Vertices)
                            {
                                if (Math.Abs(thisVaultPlanes[i].DistanceTo(n.Location)) < 1e-3)
                                {
                                    // Vertex lies directly on the plane
                                    pointsA.Add(n.Location);
                                    //Debug.WriteLine($"z: " + Math.Abs(tPlanes[i].Value.DistanceTo(v.Location)));
                                }
                                else
                                    pointsB.Add(n.Location);
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

                            foreach (var n in brep.Vertices)
                            {
                                if (Math.Abs(thisVaultPlanes[i].DistanceTo(n.Location)) < 0.001)
                                    pointsA.Add(n.Location);
                                else
                                    pointsB.Add(n.Location);
                            }
                            //Debug.WriteLine($"pointsA: " + pointsA.Count);
                            extradospoints1.AddRange(pointsA);
                            extradospoints2.AddRange(pointsB);
                        }

                        for (int j = 0; j < 2; j++) // two sides per voussoir row
                        {
                            // --- Compute new pt1 as intersection --- 
                            List<Point3d> extrados1 = (j == 0) ? extradospoints1 : extradospoints2;

                            // Find highest & lowest points
                            Point3d high = Point3d.Unset;
                            Point3d low = Point3d.Unset;
                            double maxhi = double.MinValue;
                            double minhi = double.MaxValue;

                            foreach (var p in extrados1)
                            {
                                if (p.Z > maxhi) { maxhi = p.Z; high = p; }
                                if (p.Z < minhi) { minhi = p.Z; low = p; }
                            }

                            // Build line (convert to LineCurve)
                            Line dropLine = new Line(high, low);
                            Vector3d dir = dropLine.Direction;
                            dir.Unitize(); // normalize

                            Point3d newStart = dropLine.From - dir * 300.0;
                            Point3d newEnd = dropLine.To + dir * 300.0;

                            // Rebuild extended line
                            dropLine = new Line(newStart, newEnd);

                            // Convert to curve
                            LineCurve dropCurve = new LineCurve(dropLine);
                            //Debug.WriteLine($"i + j: " + i + "," + j);
                            // Springer line as LineCurve                                                                                
                            
                            Line springerLine = springerLines[i + j];
                            LineCurve springer = new LineCurve(springerLine);

                            // Intersection
                            var events = Rhino.Geometry.Intersect.Intersection.CurveCurve(
                                springer,
                                dropCurve,
                                0.001,
                                0.001
                            );

                            Point3d pt1;
                            if (events != null && events.Count > 0)
                            {
                                pt1 = events[0].PointA;
                            }
                            else
                            {
                                // fallback
                                pt1 = springerLine.PointAt(1.0);
                            }
                            auxLinePoints.Add(pt1);
                        }
                    }
                    var ptf = auxLinePoints.First();
                    var ptl = auxLinePoints.Last();
                    var auxLine = new Line(ptf, ptl);

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
                            foreach (var n in brep.Vertices)
                            {
                                if (Math.Abs(thisVaultPlanes[i].DistanceTo(n.Location)) < 1e-3)
                                {
                                    // Vertex lies directly on the plane
                                    pointsA.Add(n.Location);
                                    //Debug.WriteLine($"z: " + Math.Abs(tPlanes[i].Value.DistanceTo(n.Location)));
                                }
                                else
                                    pointsB.Add(n.Location);
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

                            foreach (var n in brep.Vertices)
                            {
                                if (Math.Abs(thisVaultPlanes[i].DistanceTo(n.Location)) < 0.001)
                                    pointsA.Add(n.Location);
                                else
                                    pointsB.Add(n.Location);
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
                            // --- Compute new pt1 as intersection --- 
                            // Convert to curve
                            LineCurve dropCurve = new LineCurve(auxLine);

                            // Springer line as LineCurve
                            Line springerLine = springerLines[i + j];
                            LineCurve springer = new LineCurve(springerLine);

                            // Intersection
                            var events = Rhino.Geometry.Intersect.Intersection.CurveCurve(
                                springer,
                                dropCurve,
                                0.001,
                                0.001
                            );

                            Point3d pt1;

                            if (events != null && events.Count > 0)
                            {
                                pt1 = events[0].PointA;
                            }
                            else
                            {
                                // fallback
                                pt1 = springerLine.PointAt(1.0);
                            }

                            // --- Bottom inner point ---
                            Point3d pt2 = springerLines[i + j].PointAt(0);

                            // --- Intrados midpoints (used to trace interior profile) ---
                            Point3d ptA = Point3d.Unset;
                            List<Point3d> ptsA = new List<Point3d>();
                            ptsA = (j == 0) ? intradospoints1.Skip(1).ToList() : intradospoints2.Skip(1).ToList();
                            Debug.WriteLine($"ptsA: " + ptsA.Count);
                            ptA = ptsA.Last();
                            

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

                            //==============================
                            // Build polyline through ordered points
                            //==============================
                            var sPoints = new List<GH_Point>();
                            sPoints.Add(new GH_Point(pt1));
                            sPoints.Add(new GH_Point(pt2));
                            sPoints.Add(new GH_Point(ptA));
                            sPoints.Add(new GH_Point(pt5));
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

                GH_Structure<GH_Brep> springersSingle = springers;
                GH_Structure<GH_Brep> restSingle = restVoussoirs;

                // --------------------------------------------
                // 4. PREFIX OUTPUT PATHS WITH VAULT INDEX
                // --------------------------------------------
                foreach (GH_Path p in springersSingle.Paths)
                {
                    int[] newIndices = (new int[] { v }).Concat(p.Indices).ToArray();
                    GH_Path newP = new GH_Path(newIndices);

                    // Converter IList para IEnumerable<GH_Brep>
                    var branch = springersSingle.get_Branch(p).Cast<GH_Brep>();

                    outSpringers.AppendRange(branch, newP);
                }

                foreach (GH_Path p in restSingle.Paths)
                {
                    // Cria um array com v no início seguido dos índices originais
                    int[] newIndices = (new int[] { v }).Concat(p.Indices).ToArray();

                    GH_Path newP = new GH_Path(newIndices);

                    // Converter IList para IEnumerable<GH_Brep> (ou outro tipo correto) se necessário
                    var branch = restSingle.get_Branch(p).Cast<GH_Brep>();

                    outRest.AppendRange(branch, newP);
                }
            }


            // --------------------------------------------
            // 5. SET OUTPUTS
            // --------------------------------------------
            DA.SetDataTree(0, outSpringers);
            DA.SetDataTree(1, outRest);
        }
    }
}
