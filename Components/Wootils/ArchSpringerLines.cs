using Eto.Forms;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
using Rhino;
using Rhino.Geometry;
using Rhino.Geometry.Intersect;
using System;
using System.Collections.Generic;
using System.Linq;

namespace VoussoirPlugin03.Components.Wootils
{
    public class ArchSpringersFromLine : GH_Component
    {
        public ArchSpringersFromLine()
          : base("Arch Springer Lines", "ArchL",
            "Creates the springer lines for building an arch from a longitudinal curve.",
            "Voussoir", "Wootils")
        {
        }

        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                return VoussoirPlugin03.Properties.Resources.ArchSpringers;
            }
        }
        public override Guid ComponentGuid => new Guid("A1B2C3D4-F5F6-47A8-9B0C-123456789AB1");

        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddCurveParameter("Longitudinal Curve", "L", "Base Curve which describes the span of the arch.", GH_ParamAccess.item);
            pManager.AddNumberParameter("Width", "W", "Width of the arch", GH_ParamAccess.item, 0.2);
        }

        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddLineParameter("SpringerLines", "L", "Resulting springer lines", GH_ParamAccess.item);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
        // Write your logic here
        Curve C = null;
        DA.GetData(0, ref C);
        // C = C.Domain(0, 1);

        double W = 0.2;
        DA.GetData(1, ref W);


        var crvPln = new Plane(C.PointAtStart, C.PointAtEnd, C.PointAtStart + Plane.WorldXY.ZAxis);

        var Pt01 = C.PointAtStart + crvPln.ZAxis * W;
        var Pt02 = C.PointAtStart + crvPln.ZAxis * - W;

        var Pt03 = C.PointAtEnd  + crvPln.ZAxis * W;
        var Pt04 = C.PointAtEnd + crvPln.ZAxis * -W;

        var Springer01 = new Line(Pt01, Pt02);
        var Springer02 = new Line(Pt03, Pt04);

        var springer = new List<Line>();

        springer.Add(Springer01);
        springer.Add(Springer02);

        DA.SetDataList(0, springer);
        }
    }
}