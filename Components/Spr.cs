using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace VoussoirPlugin03
{
    using Components; // Ensure this is present to access Vault
    using Grasshopper;
    using Grasshopper.Kernel;
    using Grasshopper.Kernel.Data;
    using Grasshopper.Kernel.Geometry.Voronoi;
    using Grasshopper.Kernel.Parameters;
    using Grasshopper.Kernel.Types;
    using Rhino;
    using Rhino.Geometry;
    using Rhino.Geometry.Intersect;
    using System;
    using System.Collections.Generic;
    using System.Linq;
    using VoussoirPlugin03.Properties;

    namespace Components
    {
        public class VoussoirCreate : GH_Component
        {
            public VoussoirCreate()
             : base(
                   "VoussoirSpringer",
                   "Spr",
                   "Creates a Springer based on the voussoirs closest to the springer line",
                   "Voussoir",
                   "4.Springers"
                   )
            { }

            public override Guid ComponentGuid => new Guid("EC88F9F2-CD3B-4C41-ADFF-FD189794137C");

            protected override System.Drawing.Bitmap Icon
            {
                get
                {
                    // Directly return the Bitmap resource, no need to use MemoryStream
                    return VoussoirPlugin03.Properties.Resources.springer;
                }
            }
            protected override void RegisterInputParams(GH_InputParamManager pManager)
            {
                pManager.AddCurveParameter("Springer Line", "SpringerLine", "List of base lines to create vault springers", GH_ParamAccess.tree);
                pManager.AddBrepParameter("Voussoirs", "Voussoirs", "Voussoirs to analyse", GH_ParamAccess.tree);
                pManager.AddPlaneParameter("Transversal Planes", "TransversalPlanes", "Planes at each span division", GH_ParamAccess.list);
                pManager.AddNumberParameter("Springer Width", "SpringerWidth", "Distance perpendicular to springer line", GH_ParamAccess.item, 0.3);

                //Set IntradosPlanes parameter to Grafted mode
                var intradosPlanesParam = pManager[0] as Param_Curve;
                if (intradosPlanesParam != null)
                {
                    intradosPlanesParam.DataMapping = GH_DataMapping.Graft;
                }
            }

            protected override void RegisterOutputParams(GH_OutputParamManager pManager)
            {
                pManager.AddBrepParameter("Springers", "S", "Finished Springers", GH_ParamAccess.list);
                pManager.AddBrepParameter("Voussoirs", "V", "Non transformed voussoirs", GH_ParamAccess.tree);
                //pManager.AddCurveParameter("Log", "L", "All messages generated during execution", GH_ParamAccess.tree);
            }

            protected override void SolveInstance(IGH_DataAccess DA)
            {
                //=========================================
                // Declare and initialize input variables
                //=========================================
                GH_Structure<GH_Curve> spLine = new GH_Structure<GH_Curve>();
                GH_Structure<GH_Brep> voussoirs = new GH_Structure<GH_Brep>();
                List<Plane> tPlanes = new List<Plane>();
                double spWidth = 0.3;

                //=========================================
                // Retrieve input data from Grasshopper
                //=========================================
                if (!DA.GetDataTree(0, out spLine)) return;
                if (!DA.GetDataTree(1, out voussoirs)) return;
                if (!DA.GetDataList(2, tPlanes)) return;
                if (!DA.GetData(3, ref spWidth)) return;

                //=========================================
                // Extract extrados faces from voussoirs
                //=========================================
                GH_Structure<GH_Brep> extrados = new GH_Structure<GH_Brep>();

                foreach (GH_Brep b in voussoirs.AllData(true))
                {
                    Brep brep = b.Value;

                    if (brep.Faces.Count > 2)
                    {
                        BrepFace extradosFace = brep.Faces[1];
                        Brep extradosBrep = extradosFace.DuplicateFace(false);
                        extrados.Append(new GH_Brep(extradosBrep), new GH_Path(0));
                    }
                }

                //=========================================
                // Find closest extrados point to springer start, project to XY, and measure width
                //=========================================
                double smallestWidth = double.MaxValue;

                foreach (GH_Curve ghCurve in spLine.AllData(true))
                {
                    Curve springerLine = ghCurve.Value;
                    Point3d startPt = springerLine.PointAtStart;

                    // Collect all extrados points
                    List<Point3d> extradosPoints = new List<Point3d>();
                    foreach (GH_Brep ghExtrados in extrados.AllData(true))
                    {
                        Brep extradosBrep = ghExtrados.Value;
                        foreach (BrepVertex v in extradosBrep.Vertices)
                        {
                            extradosPoints.Add(v.Location);
                        }
                    }

                    // Find closest extrados point
                    Point3d closestPt = extradosPoints.OrderBy(pt => pt.DistanceTo(startPt)).FirstOrDefault();

                    // Project closest point onto XY plane
                    Point3d projectedPt = new Point3d(closestPt.X, closestPt.Y, 0);

                    // Create line and measure distance
                    double width = startPt.DistanceTo(projectedPt);
                    if (width < smallestWidth)
                        smallestWidth = width;
                }

                //=========================================
                // Create wall geometry from springer lines
                //=========================================
                List<Brep> walls = new List<Brep>();

                foreach (GH_Curve ghCurve in spLine.AllData(true))
                {
                    Curve l = ghCurve.Value;
                    Vector3d zMove = Vector3d.ZAxis * 10.0;
                    Point3d wallTopStart = l.PointAtStart + zMove;
                    Point3d wallTopEnd = l.PointAtEnd + zMove;
                    Line wallTopLine = new Line(wallTopStart, wallTopEnd);
                    Curve wallTop = wallTopLine.ToNurbsCurve();
                    Brep[] wallBreps = Brep.CreateFromLoft(new Curve[] { l, wallTop }, Point3d.Unset, Point3d.Unset, LoftType.Normal, false);
                    if (wallBreps != null && wallBreps.Length > 0)
                        walls.Add(wallBreps[0]);
                }

                //=========================================
                // Intersect walls with voussoirs and collect intersected voussoirs
                //=========================================
                GH_Structure<GH_Brep> intersectedVoussoirs = new GH_Structure<GH_Brep>();

                foreach (var path in voussoirs.Paths)
                {
                    var branch = voussoirs.get_Branch(path);
                    foreach (GH_Brep ghVoussoir in branch)
                    {
                        Brep voussoir = ghVoussoir.Value;
                        bool isIntersected = false;
                        foreach (Brep wall in walls)
                        {
                            var intersection = Rhino.Geometry.Intersect.Intersection.BrepBrep(wall, voussoir, 0.01, out var curves, out var points);
                            if (intersection && ((curves != null && curves.Length > 0) || (points != null && points.Length > 0)))
                            {
                                isIntersected = true;
                                break;
                            }
                        }
                        if (isIntersected)
                        {
                            intersectedVoussoirs.Append(ghVoussoir, path);
                        }
                    }
                }



                //=========================================
                // Output wall breps and intersected voussoirs
                //=========================================
                DA.SetDataList(0, walls);
                DA.SetDataTree(1, intersectedVoussoirs);

            }
        }
    }
}
