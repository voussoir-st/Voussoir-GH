using GH_IO.Types;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
using Rhino;
using Rhino.Geometry;
using Rhino.Geometry.Intersect;
using System;
using System.Collections.Generic;
using System.Linq;

namespace VoussoirPlugin03.Components.exp01
{
    public class centroid : GH_Component
    {
        public centroid()
         : base(
               "Geometry Centroid",
               "GCen",
               "Outputs Geometry Centroid.",
               "Voussoir",
               "Wootils"
               )
        { }

        public override Guid ComponentGuid => new Guid("F288D9F2-CD3B-4C41-ADFF-FD189794137C");

        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                // Directly return the Bitmap resource, no need to use MemoryStream
                return VoussoirPlugin03.Properties.Resources.cube09;
            }
        }
        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddGeometryParameter("Base Geometry", "G", "Geometry to calculate center", GH_ParamAccess.list);
        }

        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddPointParameter ("Centroid", "C", "The point in the bounding box center", GH_ParamAccess.list);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            List<GeometryBase> basegeometry = new List<GeometryBase>();
            List<Point3d> listadoponto = new List<Point3d>();
            
            DA.GetDataList(0, basegeometry);

            foreach (GeometryBase g in basegeometry)
            {
                Point3d centroidPt = g.GetBoundingBox(true).Center;
                listadoponto.Add(centroidPt);
            }

            DA.SetDataList(0, listadoponto);
        }

    }
}