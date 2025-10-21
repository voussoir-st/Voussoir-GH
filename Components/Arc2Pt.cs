using Grasshopper.Kernel;
using Rhino.Geometry;
using System;

namespace VoussoirPlugin03.Components
{
    public class CatenaryFromPoints : GH_Component
    {
        public CatenaryFromPoints()
          : base("Arc by Height", "ArcH",
            "Create an upward arc between two points and a given height (max radius constrained by point distance)",
            "Voussoir", "Wootils")
        {
        }

        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                return VoussoirPlugin03.Properties.Resources.Arc2Pt;
            }
        }
        public override Guid ComponentGuid => new Guid("A1B2C3D4-F5F6-47A8-9B0C-123456789ABC");

        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddPointParameter("Start", "A", "Start point of the arc.", GH_ParamAccess.item);
            pManager.AddPointParameter("End", "B", "End point of the arc.", GH_ParamAccess.item);
            pManager.AddNumberParameter("Height", "H", "Height of the arc.", GH_ParamAccess.item, 2.0);
        }

        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddCurveParameter("Arc", "C", "Arc curve.", GH_ParamAccess.item);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            Point3d a = Point3d.Unset;
            Point3d b = Point3d.Unset;
            double height = 2.0;

            if (!DA.GetData(0, ref a)) return;
            if (!DA.GetData(1, ref b)) return;
            if (!DA.GetData(2, ref height)) return;

            double dist = a.DistanceTo(b);
            if (height > dist / 2)
            {
                height = dist / 2;
            }

            Point3d mid = 0.5 * (a + b) + Vector3d.ZAxis * height;
            var arc = new Arc(a, mid, b);
            if (arc.IsValid)
            {
                DA.SetData(0, arc.ToNurbsCurve());
            }
        }
    }
}