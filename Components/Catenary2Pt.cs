using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;

namespace Components
{
    public class CatenaryFromPoints : GH_Component
    {
        public CatenaryFromPoints()
          : base("Catenary From Points", "Catenary2Pt",
            "Creates a catenary curve between two points with a specified height.",
            "Voussoir", "Utils")
        {
        }

        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                return VoussoirPlugin03.Properties.Resources.catenary;
            }
        }
        public override Guid ComponentGuid => new Guid("A1B2C3D4-E5F6-47A8-9B0C-123456789ABC");

        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddPointParameter("Start", "A", "Start point of the catenary.", GH_ParamAccess.item);
            pManager.AddPointParameter("End", "B", "End point of the catenary.", GH_ParamAccess.item);
            pManager.AddNumberParameter("Height", "H", "Rise (height) of the catenary.", GH_ParamAccess.item, 2.0);
        }

        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddCurveParameter("Catenary", "C", "Catenary curve.", GH_ParamAccess.item);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            Point3d start = Point3d.Unset;
            Point3d end = Point3d.Unset;
            double height = 2.0;

            if (!DA.GetData(0, ref start)) return;
            if (!DA.GetData(1, ref end)) return;
            if (!DA.GetData(2, ref height)) return;

            Curve catenary = BuildCatenary(start, end, height, 50);
            DA.SetData(0, catenary);
        }

        // You can copy your BuildCatenary method here, or make it static in a shared utility class.
        private Curve BuildCatenary(Point3d start, Point3d end, double rise, int samples)
        {
            if (samples < 3) samples = 3;
            var startXY = new Point2d(start.X, start.Y);
            var endXY = new Point2d(end.X, end.Y);
            Vector2d dirXY = endXY - startXY;
            double L = dirXY.Length;
            if (L < Rhino.RhinoMath.ZeroTolerance)
                return new LineCurve(start, end);
            dirXY.Unitize();
            var midXY = 0.5 * (startXY + endXY);
            double a = SolveCatenaryA(L, rise);
            if (double.IsNaN(a) || a <= 0)
            {
                Point3d apex = new Point3d(midXY.X, midXY.Y, Math.Max(start.Z, end.Z) + rise);
                return NurbsCurve.CreateInterpolatedCurve(new[] { start, apex, end }, 2);
            }
            var pts = new List<Point3d>(samples);
            for (int i = 0; i < samples; i++)
            {
                double t = i / (double)(samples - 1);
                double x = (t - 0.5) * L;
                double z = rise - (a * Math.Cosh(x / a) - a);
                double px = midXY.X + x * dirXY.X;
                double py = midXY.Y + x * dirXY.Y;
                double baseZ = start.Z + t * (end.Z - start.Z);
                pts.Add(new Point3d(px, py, baseZ + z));
            }
            if (pts.Count >= 2)
            {
                pts[0] = start;
                pts[pts.Count - 1] = end;
            }
            return Curve.CreateInterpolatedCurve(pts, 3);
        }

        private double SolveCatenaryA(double span, double rise)
        {
            if (rise <= 0 || span <= Rhino.RhinoMath.ZeroTolerance) return double.NaN;
            double a = Math.Max(span / 4.0, 1e-4);
            const double tol = 1e-10;
            const int maxIt = 100;
            for (int i = 0; i < maxIt; i++)
            {
                double x = (span * 0.5) / a;
                double cosh = Math.Cosh(x);
                double sinh = Math.Sinh(x);
                double f = a * cosh - a - rise;
                double df = cosh - 1.0 - (span * 0.5) * sinh / a;
                if (Math.Abs(df) < 1e-14) { a *= 1.2; continue; }
                double aNext = a - f / df;
                if (Math.Abs(aNext - a) < tol) return aNext;
                a = Math.Max(aNext, 1e-6);
            }
            return a;
        }
    }
}