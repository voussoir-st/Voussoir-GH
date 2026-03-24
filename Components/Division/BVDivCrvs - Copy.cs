//using Grasshopper.Kernel;
//using Grasshopper.Kernel.Parameters;
//using Rhino;
//using Rhino.Geometry;
//using System;
//using System.Collections.Generic;
//using VoussoirPlugin03.Properties;
//using Components ; // Ensure this is present to access Vault

//namespace Components
//{
//    public static class PolylineUtils
//    {
//        /// <summary>
//        /// Checks if a polyline is concave.
//        /// Returns true if concave, false if convex.
//        /// </summary>
//        public static bool IsConcave(Polyline poly)
//        {
//            if (!poly.IsClosed || poly.Count < 4)
//                throw new System.ArgumentException("Polyline must be closed and have at least 3 vertices.");

//            // Find polygon normal using first three vertices
//            Vector3d normal = Vector3d.CrossProduct(
//                poly[1] - poly[0],
//                poly[2] - poly[0]
//            );
//            normal.Unitize();

//            bool hasPositive = false;
//            bool hasNegative = false;

//            for (int i = 0; i < poly.Count - 1; i++) // -1 because last = first
//            {
//                Point3d prev = poly[(i - 1 + poly.Count - 1) % (poly.Count - 1)];
//                Point3d curr = poly[i];
//                Point3d next = poly[(i + 1) % (poly.Count - 1)];

//                Vector3d v1 = curr - prev;
//                Vector3d v2 = next - curr;

//                Vector3d cross = Vector3d.CrossProduct(v1, v2);
//                double dot = Vector3d.Multiply(cross, normal);

//                if (dot > 1e-8) hasPositive = true;
//                else if (dot < -1e-8) hasNegative = true;

//                // If both positive and negative cross products exist → concave
//                if (hasPositive && hasNegative)
//                    return true;
//            }

//            return false; // only convex turns
//        }
//        /// <summary>
//        /// Ensures that a surface's normal points in the same direction as World Z.
//        /// </summary>
//        public static Surface AlignNormalToWorldZ(Surface srf)
//        {
//            if (srf == null) return null;

//            // Pick a test point in the middle of the surface
//            double u = (srf.Domain(0).Min + srf.Domain(0).Max) * 0.5;
//            double v = (srf.Domain(1).Min + srf.Domain(1).Max) * 0.5;

//            Vector3d normal = srf.NormalAt(u, v);
//            normal.Unitize();

//            // World Z
//            Vector3d worldZ = Vector3d.ZAxis;

//            // Check dot product
//            if (normal * worldZ < 0)
//            {
//                // Flip if pointing opposite
//                srf = srf.Reverse(0); // 2 means reverse both U and V directions
//            }

//            return srf;
//        }
//    }
//    public class BarrelVaultDefiningCurves : GH_Component
//    {
//        /// <summary>
//        /// Each implementation of GH_Component must provide a public 
//        /// constructor without any arguments.
//        /// Category represents the Tab in which the component will appear, 
//        /// Subcategory the panel. If you use non-existing tab or panel names, 
//        /// new tabs/panels will automatically be created.
//        /// </summary>
//        public BarrelVaultDefiningCurves()
//          : base("Barrel Vault Defining Curves", "BVDefCrvs",
//            "Creates the defining curves for the creation of a barrel vault",
//            "Voussoir", "1.Vault Type")
//        {
//        }

//        /// <summary>
//        /// Registers all the input parameters for this component.
//        /// </summary>
//        public override Guid ComponentGuid => new Guid("6335BBF4-4F06-4327-89F9-FFDE1C891A79");

//        protected override System.Drawing.Bitmap Icon
//        {
//            get
//            {
//                // Use the imported Resources class directly
//                return Resources.BVDC;
//            }
//        }

//        protected override void RegisterInputParams(GH_InputParamManager pManager)
//        {
//            pManager.AddCurveParameter("Input Quad", "BasePoyline", "Closed Polilyne with 4 vertices.", GH_ParamAccess.list);
//            pManager.AddNumberParameter("Vault Height", "VaultHeight", "Arc's Height.", GH_ParamAccess.item, 2.0);
//            pManager.AddBooleanParameter("Span Direction", "SpanDirection", "0 = arcs on sides 0-1 and 2-3; 1 = arcs on sides 1-2 and 3-0.", GH_ParamAccess.item, true);
//            pManager.AddIntegerParameter("Arc Selection", "VaultProfile",
//                "Right-click to select type of curve\n\nParabola = 0\nArc = 1\nCatenary = 2\n\n... or input one of the above integers", GH_ParamAccess.item, 2);

//            var param = (Param_Integer)pManager[3];
//            param.AddNamedValue("Parabola", 0);
//            param.AddNamedValue("Arc", 1);
//            param.AddNamedValue("Catenary", 2);

//            pManager[1].Optional = true; // Vault Height
//            pManager[2].Optional = true; // Span Direction
//        }

//        protected override void RegisterOutputParams(GH_OutputParamManager p)
//        {
//            p.AddSurfaceParameter("Vault Surface", "VS", "The lofted base surface between the two arcs.", GH_ParamAccess.item);
//            p.AddCurveParameter("Vault Arcs", "VA", "The 2 generated arcs.", GH_ParamAccess.list);
//            p.AddCurveParameter("Lines", "L", "The 2 Horizontal lines (remaining sides).", GH_ParamAccess.list);
//        }

//        protected override void SolveInstance(IGH_DataAccess DA)
//        {
//            // Declare input variables
//            List<Curve> springerLines = new List<Curve>();
//            double height = 2.0;
//            bool spanDir = true;
//            int mode = 0;

//            // Retrieve input data from Grasshopper
//            if (!DA.GetDataList(0, springerLines))
//            {
//                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "At least 2 Curves needed.");
//                return;
//            }
//            DA.GetData(1, ref height);
//            DA.GetData(2, ref spanDir);
//            DA.GetData(3, ref mode);

//            List<Point3d> quadVertexes = new List<Point3d>();
//            foreach (Curve crv in springerLines)
//            {
//                quadVertexes.Add(crv.PointAtEnd);
//                quadVertexes.Add(crv.PointAtStart);
//            }

//            Circle circleFit;
//            Circle.TryFitCircleToPoints(quadVertexes, out circleFit);

//            // Sort quadVertexes along the fitted circle
//            var pointParams = new List<(Point3d pt, double t)>();
//            foreach (var pt in quadVertexes)
//            {
//                double t;
//                if (circleFit.ToNurbsCurve().ClosestPoint(pt, out t))
//                    pointParams.Add((pt, t));
//            }

//            // Sort by parameter value
//            pointParams.Sort((a, b) => a.t.CompareTo(b.t));

//            // Extract sorted points
//            List<Point3d> sortedQuadVertexes = new List<Point3d>();
//            foreach (var pair in pointParams)
//                sortedQuadVertexes.Add(pair.pt);

//            // Create a closed polyline from sortedQuadVertexes
//            Polyline quadCrv = new Polyline(sortedQuadVertexes);

//            // Ensure the polyline is closed
//            if (!quadCrv.IsClosed && quadCrv.Count > 0)
//                quadCrv.Add(quadCrv[0]);



//            if (Components.PolylineUtils.IsConcave(quadCrv))
//            {
//                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Input polyline cannot be concave.");
//                return;
//            }

//            // Validate that the input curve is closed
//            if (quadCrv == null || !quadCrv.IsClosed)
//            {
//                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "The curve must be closed.");
//                return;
//            }

//            // Try to extract the four corners of the quadrilateral
//            if (!VoussoirPlugin02.Components.Utils.TryGetQuadCorners(quadCrv.ToNurbsCurve(), out var p0, out var p1, out var p2, out var p3))
//            {
//                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Could not get four corners.");
//                return;
//            }

//            // Define the four edges of the quadrilateral
//            var edges = new (Point3d A, Point3d B)[]
//            {
//                (p0, p1), // 0
//                (p1, p2), // 1
//                (p2, p3), // 2
//                (p3, p0)  // 3
//            };

//            // Define which edges will be used for arcs and which for lines, based on span direction
//            int[][] arcPairs = new int[][]
//            {
//                new[] { 0, 2 }, // spanDir = 0 → arcs on sides 0-1 and 2-3
//                new[] { 1, 3 }  // spanDir = 1 → arcs on sides 1-2 and 3-0
//            };

//            int sd = spanDir ? 1 : 0;
//            var arcIdx = arcPairs[sd];
//            var lineIdx = sd == 0 ? new[] { 1, 3 } : new[] { 0, 2 };

//            // Clamp height for 3-point arc mode so both arcs have the same height and do not exceed the shortest span
//            if (mode == 1)
//            {
//                double dist0 = edges[arcIdx[0]].A.DistanceTo(edges[arcIdx[0]].B);
//                double dist1 = edges[arcIdx[1]].A.DistanceTo(edges[arcIdx[1]].B);
//                double minDist = Math.Min(dist0, dist1);
//                if (height > minDist / 2)
//                    height = minDist / 2;
//            }

//            // Build the two arc curves for the vault
//            var arcs = new List<Curve>(2);
//            foreach (int ei in arcIdx)
//            {
//                var (A, B) = edges[ei];
//                var c = BuildSpanCurve(A, B, height, mode);
//                if (c != null) arcs.Add(c);
//            }

//            // Build the two horizontal line curves for the vault
//            var lines = new List<Curve>(2);
//            foreach (int ei in lineIdx)
//            {
//                var (A, B) = edges[ei];
//                lines.Add(new LineCurve(A, B));
//            }

//            // Orient the arcs using the utility function (may flip direction for consistency)
//            VoussoirPlugin02.Components.Utils.OrientArcs(arcs);

//            // Create the lofted surface between the two arcs
//            var loftBreps = Brep.CreateFromLoft(arcs, Point3d.Unset, Point3d.Unset, LoftType.Normal, false);
//            if (loftBreps != null && loftBreps.Length > 0)
//            {
//                // Align the normal of the first lofted surface to World Z
//                var alignedSurface = Components.PolylineUtils.AlignNormalToWorldZ(loftBreps[0].Surfaces[0]);
//                Surface loftSurface = alignedSurface;
//                DA.SetData(0, loftSurface); // Output index 2 for the surface
//            }
//            else
//            {
//                DA.SetData(0, null);
//            }

//            // Output the arcs and lines to Grasshopper
//            DA.SetDataList(1, arcs); // VA: only arcs
//            DA.SetDataList(2, lines); // L: only lines
//        }

//        private Curve BuildSpanCurve(Point3d start, Point3d end, double h, int mode)
//        {
//            if (h <= RhinoMath.ZeroTolerance)
//                return new LineCurve(start, end); // altura nula → linha

//            switch (mode)
//            {
//                case 1:
//                    // Arco 3 pontos
//                    {
//                        double dist = start.DistanceTo(end);
//                        if (h > dist)
//                            h = dist; // Clamp height to the arc's span

//                        Point3d mid = 0.5 * (start + end) + Vector3d.ZAxis * h;
//                        var arc = new Arc(start, mid, end);
//                        if (arc.IsValid) return arc.ToNurbsCurve();
//                        return NurbsCurve.CreateInterpolatedCurve(new[] { start, mid, end }, 2);
//                    }
//                case 2:
//                    // Catenária com extremos a 0 e meio a +h
//                    return BuildCatenary(start, end, h, 50);
//                case 0:
//                default:
//                    // Parábola por 3 pontos
//                    {
//                        Point3d mid = 0.5 * (start + end) + Vector3d.ZAxis * h;
//                        return NurbsCurve.CreateInterpolatedCurve(new[] { start, mid, end }, 2);
//                    }
//            }
//        }

//        private Curve BuildCatenary(Point3d start, Point3d end, double rise, int samples)
//        {
//            // Draw catenary in a vertical plane between start and end (Z axis is up)
//            if (samples < 3) samples = 3;

//            // Project start and end to XY, but keep their Z values for later
//            var startXY = new Point2d(start.X, start.Y);
//            var endXY = new Point2d(end.X, end.Y);
//            Vector2d dirXY = endXY - startXY;
//            double L = dirXY.Length;
//            if (L < RhinoMath.ZeroTolerance)
//                return new LineCurve(start, end);
//            dirXY.Unitize();

//            // Find midpoint in XY
//            var midXY = 0.5 * (startXY + endXY);

//            // Solve catenary parameter 'a'
//            double a = SolveCatenaryA(L, rise);
//            if (double.IsNaN(a) || a <= 0)
//            {
//                // fallback: parabola in vertical plane
//                Point3d apex = new Point3d(midXY.X, midXY.Y, Math.Max(start.Z, end.Z) + rise);
//                return NurbsCurve.CreateInterpolatedCurve(new[] { start, apex, end }, 2);
//            }

//            var pts = new List<Point3d>(samples);
//            for (int i = 0; i < samples; i++)
//            {
//                double t = i / (double)(samples - 1); // 0..1
//                double x = (t - 0.5) * L;             // -L/2 .. L/2
//                double z = rise - (a * Math.Cosh(x / a) - a); // vertical displacement
//                // Interpolate XY position
//                double px = midXY.X + x * dirXY.X;
//                double py = midXY.Y + x * dirXY.Y;
//                // Interpolate Z between start and end, then add catenary rise
//                double baseZ = start.Z + t * (end.Z - start.Z);
//                pts.Add(new Point3d(px, py, baseZ + z));
//            }
//            // Set exact start/end points
//            if (pts.Count >= 2)
//            {
//                pts[0] = start;
//                pts[pts.Count - 1] = end;
//            }
//            return Curve.CreateInterpolatedCurve(pts, 3);
//        }

//        private double SolveCatenaryA(double span, double rise)
//        {
//            // Resolve f(a) = a*cosh((span/2)/a) - a - rise = 0
//            if (rise <= 0 || span <= RhinoMath.ZeroTolerance) return double.NaN;
//            double a = Math.Max(span / 4.0, 1e-4); // palpite inicial conservador
//            const double tol = 1e-10;
//            const int maxIt = 100;
//            for (int i = 0; i < maxIt; i++)
//            {
//                double x = (span * 0.5) / a;
//                double cosh = Math.Cosh(x);
//                double sinh = Math.Sinh(x);
//                double f = a * cosh - a - rise;
//                double df = cosh - 1.0 - (span * 0.5) * sinh / a; // derivada
//                if (Math.Abs(df) < 1e-14) { a *= 1.2; continue; }
//                double aNext = a - f / df;
//                if (Math.Abs(aNext - a) < tol) return aNext;
//                a = Math.Max(aNext, 1e-6);
//            }
//            return a; // devolve melhor estimativa
//        }
//    }
//}

////=========================================
///=========================================
//// Springer bits
///=========================================
////=========================================


////=========================================
//// 4. Create volumes from springer surfaces and Z-lines
////=========================================
//GH_Structure<GH_Brep> springerVolumes = new GH_Structure<GH_Brep>();

//for (int s = 0; s < springerSurfacesLocal.Paths.Count; s++)
//{
//    var surf = springerSurfacesLocal.Branches[s][0].Value;
//    var hLine = secondZlineLocal.Branches[s][0].Value;

//    // Compute height vector
//    Vector3d moveVec = Vector3d.ZAxis * hLine.Length;

//    // Duplicate and move surface
//    Brep surf2 = surf.DuplicateBrep();
//    surf2.Translate(moveVec);

//    // Create lofted walls between original and moved surface edges
//    var openBox = surf.Edges
//        .Select((e, i) => Brep.CreateFromLoft(
//            new List<Curve> { e.EdgeCurve, surf2.Edges[i].EdgeCurve },
//            Point3d.Unset, Point3d.Unset, LoftType.Normal, false)?.FirstOrDefault())
//        .Where(b => b != null)
//        .ToList();

//    // Combine all parts: top, bottom, and walls
//    var box = new List<Brep> { surf, surf2 };
//    box.AddRange(openBox);

//    // Join all into a single Brep volume
//    Brep joinedBox = Brep.JoinBreps(box, RhinoMath.ZeroTolerance)?.FirstOrDefault();
//    if (joinedBox != null)
//        springerVolumes.Append(new GH_Brep(joinedBox), new GH_Path(s));
//        almostSpringers.Append(new GH_Brep(joinedBox), new GH_Path(s));
//}

////=========================================
//// 5. Boolean union between volumes and intersected voussoirs
////=========================================
//for (int i = 0; i < springerVolumes.Paths.Count; i++)
//{
//    var springerBranch = springerVolumes.Branches[i];
//    var voussoirBranch = i < intersectedVoussoirs.Branches.Count ? intersectedVoussoirs.Branches[i] : null;

//    // Collect all Breps to union
//    var toUnion = springerBranch.Select(b => b?.Value)
//                    .Where(b => b != null)
//                    .ToList();

//    if (voussoirBranch != null)
//        toUnion.AddRange(voussoirBranch.Select(b => b?.Value).Where(b => b != null));

//    if (toUnion.Count == 0) continue;

//    try
//    {
//        // Attempt Boolean union
//        var unionResult = Brep.CreateBooleanUnion(toUnion, RhinoMath.ZeroTolerance);
//        if (unionResult != null && unionResult.Length > 0)
//        {
//            foreach (var b in unionResult)
//            {
//                b.MergeCoplanarFaces(RhinoMath.ZeroTolerance);
//                almostSpringers.Append(new GH_Brep(b), springerVolumes.Paths[i]);
//            }
//        }
//        else
//        {
//            // Fallback: just merge coplanar faces of original Breps
//            foreach (var b in toUnion)
//            {
//                b.MergeCoplanarFaces(RhinoMath.ZeroTolerance);
//                almostSpringers.Append(new GH_Brep(b), springerVolumes.Paths[i]);
//            }
//        }
//    }
//    catch (Exception ex)
//    {
//        RhinoApp.WriteLine($"Boolean union failed at path {springerVolumes.Paths[i]}: {ex.Message}");
//        foreach (var b in toUnion)
//        {
//            b.MergeCoplanarFaces(RhinoMath.ZeroTolerance);
//            almostSpringers.Append(new GH_Brep(b), springerVolumes.Paths[i]);
//        }
//    }
//}
//}
////==============================
//// Build Voussoirs
////==============================

//GH_Structure<GH_Brep> matchingFacesTree = new GH_Structure<GH_Brep>();

//for (int i = 0; i < almostSpringers.PathCount; i++)
//{
//    var path = almostSpringers.get_Path(i);
//    var branch = almostSpringers.Branches[i];

//    foreach (GH_Brep g in almostSpringers.Branches[i])
//    {
//        Brep brep = g.Value; // Get the actual Brep

//        // Filter faces whose planes are contained in tPlanes
//        var matchingFaces = brep.Faces
//            .Cast<BrepFace>()
//            .Where(f =>
//            {
//                Plane facePlane;
//                if (f.TryGetPlane(out facePlane))
//                {
//                    // check against all planes in tPlanes
//                    foreach (var ghPlane in tPlanes) // tPlanes is List<GH_Plane>
//                    {
//                        Plane p = ghPlane.Value; // Extract the Rhino.Geometry.Plane from GH_Plane
//                        // Normals are parallel
//                        if (facePlane.Normal.IsParallelTo(p.Normal, RhinoMath.ZeroTolerance) != 0)
//                        {
//                            // Distance from origin of one plane to the other plane is small
//                            if (Math.Abs(facePlane.DistanceTo(p.Origin)) < Rhino.RhinoMath.ZeroTolerance)
//                            {
//                                return true;
//                            }
//                        }
//                    }
//                }
//                return false;
//            })
//            .Select(f => new GH_Brep(f.DuplicateFace(false)))
//            .ToList();

//        foreach (var face in matchingFaces)
//        {
//            matchingFacesTree.Append(face, path);
//        }
//    }
//}

//int partitionSize = 2;
//GH_Structure<GH_Brep> partitionedfacesTree = TreeUtils.PartitionTree(matchingFacesTree, partitionSize);
//partitionedfacesTree.Graft(GH_GraftMode.GraftAll);

//GH_Structure<GH_Curve> edgesTree = new GH_Structure<GH_Curve>();

//Debug.WriteLine($"Original tree has {edgesTree.PathCount} paths");


//GH_Structure <GH_Point> edgesPointsTree = new GH_Structure<GH_Point>();

//for (int i = 0; i < partitionedfacesTree.PathCount; i++)
//{
//    var path = partitionedfacesTree.get_Path(i);
//    var branch = partitionedfacesTree.Branches[i];

//    foreach (var ghSurface in branch)
//    {
//        Brep brep = ghSurface.Value;
//        Curve[] edgeCurves = brep.DuplicateEdgeCurves(false);

//        var joinedCurves = Curve.JoinCurves(edgeCurves);
//        foreach (var curve in joinedCurves)
//        {
//            edgesTree.Append(new GH_Curve(curve), path);
//        }
//        edgesPointsTree.AppendRange(
//            brep.Vertices.Select(v => new GH_Point(v.Location)),
//            path
//        );
//    }
//}

//GH_Structure<GH_Curve> trimmedEdgesTree = TreeUtils.TrimTreeDepth(edgesTree);
//Debug.WriteLine($"Trimmed tree has {trimmedEdgesTree.PathCount} paths");

//var dupSpLines = new GH_Structure<GH_Curve>();

//for (int i = 0; i < matchingFacesTree.PathCount; i++)
//{
//    dupSpLines.AppendRange(spLines.Select(c => new GH_Curve(c)));
//}

//GH_Structure<GH_Curve> dupSpLinesPTree = TreeUtils.PartitionTree(dupSpLines, 2);
//dupSpLinesPTree.Simplify(GH_SimplificationMode.CollapseAllOverlaps);
//dupSpLinesPTree.Graft(GH_GraftMode.GraftAll);

//GH_Structure<GH_Curve> dupSpLinesPTree2 = TreeUtils.DuplicateBranchElements(dupSpLinesPTree);
//dupSpLinesPTree2.Graft(GH_GraftMode.GraftAll);

//GH_Structure<GH_Point> closestPointsTree = new GH_Structure<GH_Point>();
//GH_Structure<GH_Number> closestParamsTree = new GH_Structure<GH_Number>();
//GH_Structure<GH_Number> distancesTree = new GH_Structure<GH_Number>();

//TreeUtils.CurveClosestPointTree(edgesPointsTree, dupSpLinesPTree2, out closestPointsTree, out closestParamsTree, out distancesTree);
//GH_Structure<GH_Point> sortedPointsTree = TreeUtils.SortListTree(distancesTree, closestPointsTree);
//GH_Structure<GH_Point> firstPointsTree = new GH_Structure<GH_Point>();

//for (int i = 0; i < sortedPointsTree.PathCount; i++)
//{
//    var path = sortedPointsTree.get_Path(i);
//    var branch = sortedPointsTree.Branches[i];

//    if (branch.Count > 0)
//        firstPointsTree.Append(branch[0], path);
//}

//GH_Structure<GH_Point> trimmed = TreeUtils.TrimTreeDepth(firstPointsTree);

//GH_Structure<GH_Point> a = new GH_Structure<GH_Point>();
//GH_Structure<GH_Number> pr = new GH_Structure<GH_Number>();
//GH_Structure<GH_Number> bc = new GH_Structure<GH_Number>();

//TreeUtils.CurveClosestPointTree(trimmed, trimmedEdgesTree, out a, out pr, out bc);

//GH_Structure<GH_Curve> alignedEdgesTree = new GH_Structure<GH_Curve>();

//for (int i = 0; i < trimmedEdgesTree.PathCount; i++)
//{
//    var path = trimmedEdgesTree.get_Path(i);
//    var curveBranch = trimmedEdgesTree.Branches[i];

//    // Skip branches that don't have corresponding parameter branches
//    if (i >= pr.PathCount) continue;
//    var paramBranch = pr.Branches[i];

//    for (int j = 0; j < curveBranch.Count; j++)
//    {
//        var ghCurve = curveBranch[j];
//        Curve crv = ghCurve.Value;
//        if (crv == null) continue;

//        double t = 0.0;
//        if (j < paramBranch.Count)
//        {
//            t = paramBranch[j].Value;
//        }
//        else
//        {
//            // fallback: use midpoint if no parameter found
//            t = crv.Domain.Mid;
//        }

//        // Adjust seam only if curve is closed
//        if (crv.IsClosed)
//        {
//            Curve dup = crv.DuplicateCurve();
//            bool success = dup.ChangeClosedCurveSeam(t);
//            if (success)
//                alignedEdgesTree.Append(new GH_Curve(dup), path);
//            else
//                alignedEdgesTree.Append(new GH_Curve(crv.DuplicateCurve()), path);
//        }
//        else
//        {
//            alignedEdgesTree.Append(new GH_Curve(crv.DuplicateCurve()), path);
//        }
//    }
//}

//inVoussoirs.Graft(GH_GraftMode.GraftAll);

//GH_Structure<GH_Point> voussoirPoints = new GH_Structure<GH_Point>();

//for (int i = 0; i < inVoussoirs.PathCount; i++)
//{
//    var path = inVoussoirs.get_Path(i);
//    var branch = inVoussoirs.Branches[i];

//    for (int j = 0; j < branch.Count; j++)
//    {
//        var voussoirbr = branch[j];
//        var voussoir = voussoirbr.Value;
//        var voussoirpts = voussoir.Vertices;
//        voussoirPoints.AppendRange(
//            voussoirpts.Select(v => new GH_Point(v.Location)),
//            path
//        );
//    }
//}

//==============================
// Separate rest voussoirs
//==============================
//var intersectedData = inVoussoirs.AllData(true).OfType<GH_Brep>()
//    .Where(b => b.Value != null)
//    .Select(b => (Volume: VolumeMassProperties.Compute(b.Value)?.Volume ?? 0, Centroid: AreaMassProperties.Compute(b.Value)?.Centroid ?? Point3d.Unset))
//    .ToList();

//for (int i = 0; i < remappedVoussoirs.Branches.Count; i++)
//{
//    var path = remappedVoussoirs.Paths[i];
//    foreach (var ghBrep in remappedVoussoirs.Branches[i])
//    {
//        if (!(ghBrep is GH_Brep brepGoo) || brepGoo.Value == null) continue;
//        double vol = VolumeMassProperties.Compute(brepGoo.Value)?.Volume ?? 0;
//        Point3d cen = AreaMassProperties.Compute(brepGoo.Value)?.Centroid ?? Point3d.Unset;

//        if (!intersectedData.Any(d => Math.Abs(d.Volume - RhinoMath.ZeroTolerance) < RhinoMath.ZeroTolerance && d.Centroid.DistanceTo(cen) < RhinoMath.ZeroTolerance))
//            restVoussoirs.Append(ghBrep, path);
//    }
//}
