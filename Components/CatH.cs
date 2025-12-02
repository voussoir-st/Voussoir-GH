using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;

namespace Components
{
    public class CatenaryFromPoints : GH_Component
    {
        public CatenaryFromPoints()
          : base("Catenary by Height", "CatH",
            "Create an upward parabola between two points and a given height.",
            "Voussoir", "Wootils")
        {
        }

        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                return VoussoirPlugin03.Properties.Resources.CatH;
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
            //pManager.AddCurveParameter("Log", "L", "All messages generated during execution", GH_ParamAccess.tree);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {


            Point3d start = Point3d.Unset;
            Point3d end = Point3d.Unset;
            double height = 2.0;

            if (!DA.GetData(0, ref start)) return;
            if (!DA.GetData(1, ref end)) return;
            if (!DA.GetData(2, ref height)) return;

            Curve catenary = BuildCatenary(start, end, height, 30);
            Curve a = catenary.Rebuild(20, 3, true);
            DA.SetData(0, catenary);
        }

        private Curve BuildCatenary(Point3d start, Point3d end, double rise, int samples)
        {
            if (samples < 3) samples = 3;
            
            Vector3d dirXY = end - start;
            double L = start.DistanceTo(end);
            if (L < Rhino.RhinoMath.ZeroTolerance)
                return new LineCurve(start, end);
            dirXY.Unitize();
            var mid = 0.5 * (start + end);
            double a = SolveCatenaryA(L, rise);
            
            var pts = new List<Point3d>(samples);
            for (int i = 0; i < samples; i++)
            {
                double t = i / (double)(samples - 1);
                double x = (t - 0.5) * L;
                double z = rise - (a * Math.Cosh(x / a) - a);
                double px = mid.X + x * dirXY.X;
                double py = mid.Y + x * dirXY.Y;
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
            //List<string> log = new List<string>();
            //log.Clear();

            double a = Math.Max(span / 10.0, 1e-14);
            //AddRuntimeMessage(GH_RuntimeMessageLevel.Remark, "debug 02 a:" + a);
            const double tol = 0;
            const int maxIt = 1000; //10000
            for (int i = 0; i < maxIt; i++)
            {
                //AddRuntimeMessage(GH_RuntimeMessageLevel.Remark, "debug 00 path:" + i);
                //log.Add("hello"+i);
                double x = (span * 0.5) / a;
                double cosh = Math.Cosh(x);
                double sinh = Math.Sinh(x);
                double f = a * cosh - a - rise;
                double df = cosh - 1.0 - (span * 0.5) * sinh / a;
                //AddRuntimeMessage(GH_RuntimeMessageLevel.Remark, "debug 00 path:" + i);
                //if (Math.Abs(df) < 1e-14) { a *= 1.2; continue; }
                //if (Math.Abs(df) < 1e-1) { a *= 1.2; continue; }
                //AddRuntimeMessage(GH_RuntimeMessageLevel.Remark, "debug 02 path:" + i);
                double aNext = a - f / df;
                if (Math.Abs(aNext - a) < tol) return aNext;
                a = Math.Max(aNext, 1e-2);
                //AddRuntimeMessage(GH_RuntimeMessageLevel.Remark, "debug 03 a:" + a);
            }
            //AddRuntimeMessage(GH_RuntimeMessageLevel.Remark, "debug 04 a:" + a);
            return a;

        }
    }
}