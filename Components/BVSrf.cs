using Components ;
using Ed.Eto; // Ensure this is present to access Vault
using Grasshopper.Kernel;
using Grasshopper.Kernel.Parameters;
using Rhino;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using VoussoirPlugin03.Properties;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
using Rhino;
using Rhino.Geometry;
using Rhino.Geometry.Intersect;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Security.Cryptography;
using VoussoirPlugin03.Properties;

namespace Components
{
    public static class PolylineUtils
    {
        public static Vector3d GetBrepNormalAtPoint(Brep brep, Point3d pt, double tol)
        {
            // Find closest point on the brep
            ComponentIndex ci;
            double u, v;
            Point3d closestPt;
            Vector3d normal;

            if (!brep.ClosestPoint(pt, out closestPt, out ci, out u, out v, tol, out normal))
                return Vector3d.Unset;

            // Ensure it's a face
            if (ci.ComponentIndexType != ComponentIndexType.BrepFace)
                return Vector3d.Unset;

            BrepFace face = brep.Faces[ci.Index];

            // Get the normal
            Vector3d faceNormal = face.NormalAt(u, v);
            faceNormal.Unitize();

            return faceNormal;
        }
        /// <summary>
        /// Checks if a polyline is concave.
        /// Returns true if concave, false if convex.
        /// </summary>
        public static bool IsConcave(Polyline poly)
        {
            if (!poly.IsClosed || poly.Count < 4)
                throw new System.ArgumentException("Polyline must be closed and have at least 3 vertices.");

            // Find polygon normal using first three vertices
            Vector3d normal = Vector3d.CrossProduct(
                poly[1] - poly[0],
                poly[2] - poly[0]
            );
            normal.Unitize();

            bool hasPositive = false;
            bool hasNegative = false;

            for (int i = 0; i < poly.Count - 1; i++) // -1 because last = first
            {
                Point3d prev = poly[(i - 1 + poly.Count - 1) % (poly.Count - 1)];
                Point3d curr = poly[i];
                Point3d next = poly[(i + 1) % (poly.Count - 1)];

                Vector3d v1 = curr - prev;
                Vector3d v2 = next - curr;

                Vector3d cross = Vector3d.CrossProduct(v1, v2);
                double dot = Vector3d.Multiply(cross, normal);

                if (dot > 1e-8) hasPositive = true;
                else if (dot < -1e-8) hasNegative = true;

                // If both positive and negative cross products exist → concave
                if (hasPositive && hasNegative)
                    return true;
            }

            return false; // only convex turns
        }
        /// <summary>
        /// Ensures that a surface's normal points in the same direction as World Z.
        /// </summary>
        public static Surface AlignNormalToWorldZ(Surface srf)
        {
            if (srf == null) return null;

            // Pick a test point in the middle of the surface
            double u = (srf.Domain(0).Min + srf.Domain(0).Max) * 0.5;
            double v = (srf.Domain(1).Min + srf.Domain(1).Max) * 0.5;

            Vector3d normal = srf.NormalAt(u, v);
            normal.Unitize();

            // World Z
            Vector3d worldZ = Vector3d.ZAxis;

            // Check dot product
            if (normal * worldZ < 0)
            {
                // Flip if pointing opposite
                srf = srf.Reverse(0); // 2 means reverse both U and V directions
            }

            return srf;
        }
    }
    public class BarrelVaultDefiningCurves : GH_Component
    {
        /// <summary>
        /// Each implementation of GH_Component must provide a public 
        /// constructor without any arguments.
        /// Category represents the Tab in which the component will appear, 
        /// Subcategory the panel. If you use non-existing tab or panel names, 
        /// new tabs/panels will automatically be created.
        /// </summary>
        public BarrelVaultDefiningCurves()
          : base(
            "Barrel Vault Base Surface", 
            "BVSrf",
            "Creates the base surface for the creation of a barrel vault",
            "Voussoir",
            "Core Geometry"
          )
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        public override Guid ComponentGuid => new Guid("6335BBF4-4F06-4327-89F9-FFDE1C891A79");

        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                // Use the imported Resources class directly
                return Resources.BVSrf;
            }
        }

        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddCurveParameter("Input Lines", "SpringerLines", "set(s) of 2 non intersecting lines", GH_ParamAccess.list);
            pManager.AddNumberParameter("Vault Height", "VaultHeight", "Arc's Height.", GH_ParamAccess.item, 2.0);
            //pManager.AddBooleanParameter("Span Direction", "SpanDirection", "0 = arcs on sides 0-1 and 2-3; 1 = arcs on sides 1-2 and 3-0.", GH_ParamAccess.item, true);
            pManager.AddIntegerParameter("Arc Selection", "VaultProfile",
                "Right-click to select type of curve\n\nArc = 0\nParabola = 1\nCatenary = 2\n\n... or input one of the above integers", GH_ParamAccess.item, 0);

            var param = (Param_Integer)pManager[2];
            param.AddNamedValue("Arc", 0);
            param.AddNamedValue("Parabola", 1);
            param.AddNamedValue("Catenary", 2);

            pManager[1].Optional = true; // Vault Height
            pManager[2].Optional = true; // Span Direction
        }

        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddSurfaceParameter("Vault Surface", "VS", "The lofted base surface between the two arcs.", GH_ParamAccess.tree);
            pManager.AddCurveParameter("Springer Lines", "SL", "The 2 Horizontal lines (remaining sides).", GH_ParamAccess.tree);
            pManager.HideParameter(1);
            pManager.AddCurveParameter("Vault Arcs", "VA", "The 2 generated arcs.", GH_ParamAccess.tree);
            pManager.HideParameter(2);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // Declare input variables
            List<Curve> springerLines = new List<Curve>();
            double height = 2.0;
            bool spanDir = false;
            int mode = 0;

            // Retrieve input data from Grasshopper
            // 1. Try to get data
            if (!DA.GetDataList(0, springerLines))
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "No curves supplied.");
                return;   // stop execution here
            }

            // 2. Check list size
            if (springerLines.Count < 2)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "At least 2 curves needed.");
                return;   // stop execution
            }

            DA.GetData(1, ref height);
            //DA.GetData(2, ref spanDir);
            DA.GetData(2, ref mode);

            GH_Path currentPath = this.Params.Input[0].VolatileData.get_Path(DA.Iteration);
            var linesTree = new GH_Structure<GH_Curve>();
            var arcsTree = new GH_Structure<GH_Curve>();

            Components.Utils.OrientArcs(springerLines);
            bool intersectLines = false;
            if (springerLines.Count >= 2)
            {
                var events = Rhino.Geometry.Intersect.Intersection.CurveCurve(springerLines[0], springerLines[1], RhinoDoc.ActiveDoc.ModelAbsoluteTolerance, 0.0);
                intersectLines = events != null && events.Count > 0;

                if (intersectLines)
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Springer Lines cannot intersect.");
                    return;
                }
            }

            // Define the four edges of the quadrilateral
            var edges = new (Point3d A, Point3d B)[]
            {
                (springerLines[0].PointAtStart, springerLines[1].PointAtStart), // 0
                (springerLines[1].PointAtStart, springerLines[1].PointAtEnd), // 1
                (springerLines[1].PointAtEnd, springerLines[0].PointAtEnd), // 2
                (springerLines[0].PointAtEnd, springerLines[0].PointAtStart), // 3
            };

            // Define which edges will be used for arcs and which for lines, based on span direction
            int[][] arcPairs = new int[][]
            {
                new[] { 0, 2 }, // spanDir = 0 → arcs on sides 0-1 and 2-3
                new[] { 1, 3 }  // spanDir = 1 → arcs on sides 1-2 and 3-0
            };

            int sd = spanDir ? 1 : 0;
            var arcIdx = arcPairs[sd];
            var lineIdx = sd == 0 ? new[] { 1, 3 } : new[] { 0, 2 };

            // Clamp height for 3-point arc mode so both arcs have the same height and do not exceed the shortest span
            if (mode == 0)
            {
                double dist0 = edges[arcIdx[0]].A.DistanceTo(edges[arcIdx[0]].B);
                double dist1 = edges[arcIdx[1]].A.DistanceTo(edges[arcIdx[1]].B);
                double minDist = Math.Min(dist0, dist1);
                if (height > minDist / 2)
                    height = minDist / 2;
            }

            // Build the two arc curves for the vault
            var arcs = new List<Curve>(2);
            foreach (int ei in arcIdx)
            {
                var (A, B) = edges[ei];
                var c = BuildSpanCurve(A, B, height, mode);
                Curve d = c.Rebuild(20, 3, true);             
                if (d != null) arcs.Add(d);
            }

            // Build the two horizontal line curves for the vault
            var lines = new List<Curve>(2);
            foreach (int ei in lineIdx)
            {
                var (A, B) = edges[ei];
                lines.Add(new LineCurve(A, B));
            }

            // Orient the arcs using the utility function (may flip direction for consistency)
            Components.Utils.OrientArcs(arcs);

            // Create the lofted surface between the two arcs
            var loftBreps = Brep.CreateFromLoft(arcs, Point3d.Unset, Point3d.Unset, LoftType.Normal, false);
            if (loftBreps != null && loftBreps.Length > 0)
            {
                // Align the normal of the first lofted surface to World Z
                var alignedSurface = Components.PolylineUtils.AlignNormalToWorldZ(loftBreps[0].Surfaces[0]);
                Surface loftSurface = alignedSurface;
                DA.SetData(0, loftSurface); // Output index 2 for the surface
            }
            else
            {
                DA.SetData(0, null);
            }

            lines.Reverse();
            Components.Utils.OrientArcs(lines);

            foreach (var l in lines)
                linesTree.Append(new GH_Curve(l), currentPath);

            foreach (var a in arcs)
                arcsTree.Append(new GH_Curve(a), currentPath);

            DA.SetDataTree(1, linesTree);
            DA.SetDataTree(2, arcsTree);
        }

        private Curve BuildSpanCurve(Point3d start, Point3d end, double h, int mode)
        {
            if (h <= RhinoMath.ZeroTolerance)
                return new LineCurve(start, end); // altura nula → linha

            switch (mode)
            {
                case 0:
                    // Arco 3 pontos
                    {
                        double dist = start.DistanceTo(end);
                        if (h > dist)
                            h = dist; // Clamp height to the arc's span

                        Point3d mid = 0.5 * (start + end) + Vector3d.ZAxis * h;
                        var arc = new Arc(start, mid, end);
                        if (arc.IsValid) return arc.ToNurbsCurve();
                        return NurbsCurve.CreateInterpolatedCurve(new[] { start, mid, end }, 2);
                    }
                case 2:
                    // Catenária com extremos a 0 e meio a +h
                    return BuildCatenary(start, end, h, 50);
                case 1:
                default:
                    // Parábola por 3 pontos
                    {
                        Point3d mid = 0.5 * (start + end) + Vector3d.ZAxis * h;
                        return NurbsCurve.CreateInterpolatedCurve(new[] { start, mid, end }, 2);
                    }
            }
        }

        private Curve BuildCatenary(Point3d start, Point3d end, double rise, int samples)
        {
            // Draw catenary in a vertical plane between start and end (Z axis is up)
            if (samples < 3) samples = 3;

            // Project start and end to XY, but keep their Z values for later
            var startXY = new Point2d(start.X, start.Y);
            var endXY = new Point2d(end.X, end.Y);
            Vector2d dirXY = endXY - startXY;
            double L = dirXY.Length;
            if (L < RhinoMath.ZeroTolerance)
                return new LineCurve(start, end);
            dirXY.Unitize();

            // Find midpoint in XY
            var midXY = 0.5 * (startXY + endXY);

            // Solve catenary parameter 'a'
            double a = SolveCatenaryA(L, rise);
            if (double.IsNaN(a) || a <= 0)
            {
                // fallback: parabola in vertical plane
                Point3d apex = new Point3d(midXY.X, midXY.Y, Math.Max(start.Z, end.Z) + rise);
                return NurbsCurve.CreateInterpolatedCurve(new[] { start, apex, end }, 2);
            }

            var pts = new List<Point3d>(samples);
            for (int i = 0; i < samples; i++)
            {
                double t = i / (double)(samples - 1); // 0..1
                double x = (t - 0.5) * L;             // -L/2 .. L/2
                double z = rise - (a * Math.Cosh(x / a) - a); // vertical displacement
                // Interpolate XY position
                double px = midXY.X + x * dirXY.X;
                double py = midXY.Y + x * dirXY.Y;
                // Interpolate Z between start and end, then add catenary rise
                double baseZ = start.Z + t * (end.Z - start.Z);
                pts.Add(new Point3d(px, py, baseZ + z));
            }
            // Set exact start/end points
            if (pts.Count >= 2)
            {
                pts[0] = start;
                pts[pts.Count - 1] = end;
            }
            return Curve.CreateInterpolatedCurve(pts, 3);
        }

        private double SolveCatenaryA(double span, double rise)
        {
            // Resolve f(a) = a*cosh((span/2)/a) - a - rise = 0
            if (rise <= 0 || span <= RhinoMath.ZeroTolerance) return double.NaN;
            double a = Math.Max(span / 4.0, 1e-4); // palpite inicial conservador
            const double tol = 1e-10;
            const int maxIt = 100;
            for (int i = 0; i < maxIt; i++)
            {
                double x = (span * 0.5) / a;
                double cosh = Math.Cosh(x);
                double sinh = Math.Sinh(x);
                double f = a * cosh - a - rise;
                double df = cosh - 1.0 - (span * 0.5) * sinh / a; // derivada
                if (Math.Abs(df) < 1e-14) { a *= 1.2; continue; }
                double aNext = a - f / df;
                if (Math.Abs(aNext - a) < tol) return aNext;
                a = Math.Max(aNext, 1e-6);
            }
            return a; // devolve melhor estimativa
        }
    }
}
