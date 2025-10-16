using Components; // Ensure this is present to access Vault
using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Attributes;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Parameters;
using Grasshopper.Kernel.Types;
using Rhino;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using System.Linq;
using VoussoirPlugin03.Properties;

namespace Components
{
    public class IrregularVaultDivisionComponent : GH_Component
    {
        public IrregularVaultDivisionComponent()
          : base(
                "Irregular Vault Division",
                "IVDiv",
                "Divides a vault defined by two arcs into spanwise and lengthwise voussoirs.",
                "Voussoir",
                "2.Division"
                )
        { }

        public override Guid ComponentGuid => new Guid("BC98D9F2-CD3B-4C41-ADFF-FD189794437C");

        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                return VoussoirPlugin03.Properties.Resources.VDiv;
            }
        }
        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddSurfaceParameter("Vault Surface", "VaultSurface", "Surface that defines the vault", GH_ParamAccess.item);
            pManager.AddIntegerParameter("Transversal Divisions", "TransversalDivisions", "Number of voussoirs in the vault's span", GH_ParamAccess.item, 12);
            pManager.AddIntegerParameter("Longitudinal divisions", "LongitudinalDivisions", "Number of voussoirs in the vault's length", GH_ParamAccess.item, 8);
        }

        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddPlaneParameter("Intrados Planes", "IP", "Intrados planar vault panels", GH_ParamAccess.list);
            pManager.AddPlaneParameter("Division planes", "DP", "Planes of each Voussoir Contact Surface", GH_ParamAccess.tree);
            pManager.AddPlaneParameter("Transversal planes", "TP", "Planes at each span division", GH_ParamAccess.list);
            pManager.AddPlaneParameter("Longitudinal planes", "LP", "Planes at each length division", GH_ParamAccess.list);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            Surface vaultSrf = null;
            int spanDiv = 12;
            int lengthDiv = 8;

            if (!DA.GetData(0, ref vaultSrf)) return;
            if (!DA.GetData(2, ref spanDiv)) return;
            if (!DA.GetData(1, ref lengthDiv)) return;

            // 1. Create Lofted Surface
            Brep loftBrep = vaultSrf.ToBrep();
            if (loftBrep == null || loftBrep.Edges.Count < 2)
                return;
            Surface loftSrf = Components.PolylineUtils.AlignNormalToWorldZ(loftBrep.Surfaces[0]);

            // 2. Extract boundary curves (arcs) from the surface
            Curve[] boundaryCurvesArray = loftBrep.DuplicateEdgeCurves(true);
            List<Curve> boundaryCurves = boundaryCurvesArray.ToList();
            if (boundaryCurves.Count < 2)
                return;

            // Assume the first and last boundary curves in the V direction are the vault arcs
            Curve arc1 = boundaryCurves[1];
            Curve arc2 = boundaryCurves[3]; // This may need adjustment depending on surface orientation

            if (!arc1.IsValid || !arc2.IsValid || spanDiv < 2 | lengthDiv < 2) return;

            // Use Utils.OrientArcs if needed
            Components.Utils.OrientArcs(new List<Curve> { arc1, arc2 });
            arc1 = boundaryCurves[1];
            arc2 = boundaryCurves[3];

            // Otherwise, keep DivideByCount for count-based division
            List<Point3d> arc1_pts = new List<Point3d>();
            List<Point3d> arc2_pts = new List<Point3d>();
            double[] params1 = arc1.DivideByCount(spanDiv, true);
            double[] params2 = arc2.DivideByCount(spanDiv, true);

            foreach (var t in params1) arc1_pts.Add(arc1.PointAt(t));
            foreach (var t in params2) arc2_pts.Add(arc2.PointAt(t));

            // 3. Build grid: grid[s][l], size [spanDiv+1][lengthDiv+1]
            List<List<Point3d>> grid = new List<List<Point3d>>();
            for (int s = 0; s <= spanDiv; s++)
            {
                var row = new List<Point3d>();
                Line spanLine = new Line(arc1_pts[s], arc2_pts[s]);
                // Use Utils.DivideLineByLength if you want length-based division
                // Otherwise, keep count-based division
                for (int l = 0; l <= lengthDiv; l++)
                {
                    double factor = (double)l / lengthDiv; // 0 to 1
                    row.Add(spanLine.PointAt(factor));
                }
                grid.Add(row);
            }

            //4. building the longitudinal planes
            List<Plane> longitudinalPlanes = new List<Plane>();
            double offsetDistance = 1.0;
            for (int s = 0; s <= spanDiv; s++)
            {
                double sumX = 0, sumY = 0, sumZ = 0;
                for (int l = 0; l <= lengthDiv; l++)
                {
                    Point3d pt = grid[s][l];
                    sumX += pt.X;
                    sumY += pt.Y;
                    sumZ += pt.Z;
                }
                double avgX = sumX / (lengthDiv + 1);
                double avgY = sumY / (lengthDiv + 1);
                double avgZ = sumZ / (lengthDiv + 1);
                Point3d avgPt = new Point3d(avgX, avgY, avgZ);

                Point3d lastPt = grid[s][lengthDiv];

                List<Point3d> movedPts = new List<Point3d>();
                for (int l = 0; l <= lengthDiv; l++)
                {
                    Point3d pt = grid[s][l];
                    double u, v;
                    if (loftSrf.ClosestPoint(pt, out u, out v))
                    {
                        Vector3d normal = loftSrf.NormalAt(u, v);
                        normal.Unitize();
                        Point3d ptMoved = pt + normal * offsetDistance;
                        movedPts.Add(ptMoved);
                    }
                    else
                    {
                        movedPts.Add(pt); // fallback
                    }
                }
                double sumMX = 0, sumMY = 0, sumMZ = 0;
                foreach (var mpt in movedPts)
                {
                    sumMX += mpt.X;
                    sumMY += mpt.Y;
                    sumMZ += mpt.Z;
                }
                Point3d avgMovedPt = new Point3d(
                    sumMX / movedPts.Count,
                    sumMY / movedPts.Count,
                    sumMZ / movedPts.Count
                );

                Plane longitudinalPlane = new Plane(avgPt, lastPt, avgMovedPt);
                longitudinalPlanes.Add(longitudinalPlane);
            }
           
            //5. building the transversal planes
            List<Plane> transversalPlanes = new List<Plane>();
            for (int l = 0; l <= lengthDiv; l++)
            {
                List<Point3d> pts = new List<Point3d>();
                List<Point3d> movedPts = new List<Point3d>();
                for (int s = 0; s <= spanDiv; s++)
                {
                    Point3d pt = grid[s][l];
                    pts.Add(pt);

                    double u, v;
                    if (loftSrf.ClosestPoint(pt, out u, out v))
                    {
                        Point3d ptMoved = pt + Vector3d.ZAxis * offsetDistance;

                        movedPts.Add(ptMoved);
                    }
                    else
                    {
                        movedPts.Add(pt); // fallback
                    }
                }

                double sumX = 0, sumY = 0, sumZ = 0;
                foreach (var pt in pts)
                {
                    sumX += pt.X;
                    sumY += pt.Y;
                    sumZ += pt.Z;
                }
                Point3d avgPt = new Point3d(sumX / pts.Count, sumY / pts.Count, sumZ / pts.Count);

                Point3d lastPt = grid[spanDiv][l];

                double sumMX = 0, sumMY = 0, sumMZ = 0;
                foreach (var mpt in movedPts)
                {
                    sumMX += mpt.X;
                    sumMY += mpt.Y;
                    sumMZ += mpt.Z;
                }
                Point3d avgMovedPt = new Point3d(sumMX / movedPts.Count, sumMY / movedPts.Count, sumMZ / movedPts.Count);

                Plane transversalPlane = new Plane(avgPt, lastPt, avgMovedPt);
                transversalPlanes.Add(transversalPlane);
            }

            GH_Structure<GH_Plane> panelPlanesTree = new GH_Structure<GH_Plane>();
            Plane referencePlane = Plane.Unset;

            for (int s = 0; s < spanDiv; s++)
            {
                GH_Path path = new GH_Path(s);
                for (int l = 0; l < lengthDiv; l++)
                {
                    var a = grid[s][l];
                    var b = grid[s + 1][l];
                    var c = grid[s + 1][l + 1];
                    var d = grid[s][l + 1];
                    var pts = new Point3d[] { a, b, c, d };

                    Plane fitPlane;
                    Plane.FitPlaneToPoints(pts, out fitPlane);

                    double avgX = (a.X + b.X + c.X + d.X) / 4.0;
                    double avgY = (a.Y + b.Y + c.Y + d.Y) / 4.0;
                    double avgZ = (a.Z + b.Z + c.Z + d.Z) / 4.0;
                    Point3d avgPt = new Point3d(avgX, avgY, avgZ);

                    // Get the normal of the lofted surface at the average point
                    double u, v;
                    Vector3d surfNormal = fitPlane.ZAxis; // fallback
                    if (loftSrf.ClosestPoint(avgPt, out u, out v))
                    {
                        surfNormal = loftSrf.NormalAt(u, v);
                        surfNormal.Unitize();
                    }

                    Vector3d panelNormal = fitPlane.ZAxis;
                    panelNormal.Unitize();

                    // Flip the plane if its Z axis is opposed to the surface normal
                    if (!panelNormal.IsZero && !surfNormal.IsZero && panelNormal * surfNormal < 0)
                    {
                        fitPlane.Flip();
                        panelNormal = fitPlane.ZAxis;
                        panelNormal.Unitize();
                    }

                    double angle = Vector3d.VectorAngle(fitPlane.XAxis, Plane.WorldXY.XAxis, fitPlane);
                    if (angle > Rhino.RhinoMath.ZeroTolerance)
                    {
                        fitPlane.Rotate(angle, panelNormal);
                    }

                    Plane panelPlane = new Plane(avgPt, fitPlane.XAxis, fitPlane.YAxis);
                    panelPlanesTree.Append(new GH_Plane(panelPlane), path);
                }
            }

            // Replace the divisionPlanesTree construction section with this to ensure proper tree organization
            GH_Structure<GH_Plane> divisionPlanesTree = new GH_Structure<GH_Plane>();

            for (int s = 0; s < spanDiv; s++)
            {
                for (int l = 0; l < lengthDiv; l++)
                {
                    // Each branch is uniquely identified by {s; l}
                    GH_Path path = new GH_Path(s, l);

                    // Collect the four planes for this cell
                    divisionPlanesTree.Append(new GH_Plane(longitudinalPlanes[s]), path);     // Start longitudinal
                    divisionPlanesTree.Append(new GH_Plane(longitudinalPlanes[s + 1]), path); // End longitudinal
                    divisionPlanesTree.Append(new GH_Plane(transversalPlanes[l]), path);      // Start transversal
                    divisionPlanesTree.Append(new GH_Plane(transversalPlanes[l + 1]), path);  // End transversal
                }
            }

            DA.SetDataTree(0, panelPlanesTree);
            DA.SetDataTree(1, divisionPlanesTree);
            DA.SetDataList(3, longitudinalPlanes);
            DA.SetDataList(2, transversalPlanes);
        }

    }
}