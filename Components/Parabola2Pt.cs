using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using VoussoirPlugin03.Properties;

namespace VoussoirPlugin03.Components
{
    public class ParabolaFromPoints : GH_Component
    {
        public ParabolaFromPoints()
          : base("Parabola by Height", "ParH",
            "Create an upward catenary between two points and a given height.",
            "Voussoir", "Wootils")
        {
        }

        public override Guid ComponentGuid => new Guid("B2C3D4E5-F6A7-48B9-9C0D-23456789ABCD");

        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                // Use the imported Resources class directly
                return Resources.Parabola2Pt;
            }
        }

        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddPointParameter("Start", "A", "Start point of the parabola.", GH_ParamAccess.item);
            pManager.AddPointParameter("End", "B", "End point of the parabola.", GH_ParamAccess.item);
            pManager.AddNumberParameter("Height", "H", "Height of the parabola.", GH_ParamAccess.item, 2.0);
        }

        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddCurveParameter("Parabola", "C", "Parabola curve.", GH_ParamAccess.item);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            Point3d start = Point3d.Unset;
            Point3d end = Point3d.Unset;
            double height = 2.0;

            if (!DA.GetData(0, ref start)) return;
            if (!DA.GetData(1, ref end)) return;
            if (!DA.GetData(2, ref height)) return;

            Curve parabola = BuildParabola(start, end, height);
            DA.SetData(0, parabola);
        }

        private Curve BuildParabola(Point3d start, Point3d end, double rise)
        {
            // Midpoint between start and end, elevated by rise in Z
            Point3d mid = 0.5 * (start + end) + Vector3d.ZAxis * rise;
            // Interpolated curve through three points
            return NurbsCurve.CreateInterpolatedCurve(new[] { start, mid, end }, 2);
        }
    }
}