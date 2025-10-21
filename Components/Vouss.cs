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
               "Create Voussoirs",
               "BVouss",
               "Creates Voussoirs based on predefined planes.",
               "Voussoir",
               "Vault Creation"
               )
        { }

        public override Guid ComponentGuid => new Guid("EC88D9F2-CD3B-4C41-ADFF-FD189794137C");

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
            pManager.AddPlaneParameter("Intrados Planes", "IntradosPlanes", "Intrados planar vault panels", GH_ParamAccess.tree);
            pManager.AddPlaneParameter("Division planes", "DivisionPlanes", "Planes of each Voussoir Contact Surface", GH_ParamAccess.tree);
            pManager.AddNumberParameter("Voussoir Thickness", "Thickness", "User defined Voussoir Thickness", GH_ParamAccess.item, 0.3);

            // Set IntradosPlanes parameter to Grafted mode
            var intradosPlanesParam = pManager[0] as Param_Plane;
            if (intradosPlanesParam != null)
            {
               intradosPlanesParam.DataMapping = GH_DataMapping.Graft;
            }
        }

        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddBrepParameter("Voussoir", "V", "Voussoir Blocks", GH_ParamAccess.tree);
            pManager.AddSurfaceParameter("Extrados", "E", "Extrados Surfaces", GH_ParamAccess.tree);
            pManager.HideParameter(1);
            pManager.AddSurfaceParameter("Intrados", "I", "Intrados Surfaces", GH_ParamAccess.tree);
            pManager.HideParameter(2);
            pManager.AddBrepParameter("Contact Faces", "CF", "Voussoir contact faces", GH_ParamAccess.tree);
            pManager.HideParameter(3);
            //pManager.AddCurveParameter("Log", "L", "All messages generated during execution", GH_ParamAccess.tree);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            
            //List<string> log = new List<string>();
            //log.Clear();
            GH_Structure<GH_Plane> intradosPlanes;
            GH_Structure<GH_Plane> divPlanes;
            double thickness = 0.3;

            //log.Add("hello");

            //RhinoApp.WriteLine("SolveInstance called");

            if (!DA.GetDataTree(0, out intradosPlanes)) return;
            if (!DA.GetDataTree(1, out divPlanes)) return;
            if (!DA.GetData(2, ref thickness)) return;

            DataTree<Line> intersectionCurve = new DataTree<Line>();
            DataTree<Point3d> intersectionPoints = new DataTree<Point3d>();

            int[,] pairs = new int[,] { { 0, 2 }, { 0, 3 }, { 1, 2 }, { 1, 3 } };

            foreach (GH_Path path in divPlanes.Paths)
            {
                //AddRuntimeMessage(GH_RuntimeMessageLevel.Remark, "debug 00 path:" + path);
                var branchGoo = divPlanes.get_Branch(path);
                if (branchGoo == null || branchGoo.Count != 4)
                {
                    RhinoApp.WriteLine($"Branch {path} doesn't have 4 elements ({branchGoo?.Count ?? 0} found).");
                    continue;
                }

                List<Plane> branchPlanes = branchGoo.Cast<GH_Plane>().Select(p => p.Value).ToList();

                Plane intradosPlane = ((GH_Plane)intradosPlanes.get_Branch(path)[0]).Value;

                HashSet<Line> uniqueLines = new HashSet<Line>();
                HashSet<Point3d> uniquePoints = new HashSet<Point3d>();

                for (int k = 0; k < pairs.GetLength(0); k++)
                {

                    //AddRuntimeMessage(GH_RuntimeMessageLevel.Remark, "debug 01 k:"+k);
                    //AddRuntimeMessage(GH_RuntimeMessageLevel.Remark, "debug 02");

                    int a = pairs[k, 0];
                    int b = pairs[k, 1];

                    Plane planeA = branchPlanes[a];
                    Plane planeB = branchPlanes[b];

                    //=========================================
                    // Intersect planes
                    //=========================================

                    if (Rhino.Geometry.Intersect.Intersection.PlanePlane(planeA, planeB, out Line line))
                    {
                        //AddRuntimeMessage(GH_RuntimeMessageLevel.Remark, "debug 03 8x? k:"+k);

                        if (uniqueLines.Add(line))
                            intersectionCurve.Add(line, path);
                        //AddRuntimeMessage(GH_RuntimeMessageLevel.Remark, "adicionou line. "+path+" k:"+ k );

                        //=========================================
                        // Intersect line with intrados plane
                        //=========================================

                        if (Rhino.Geometry.Intersect.Intersection.LinePlane(line, intradosPlane, out double t))
                        {
                            Point3d point = line.PointAt(t);
                            if (uniquePoints.Add(point))
                                intersectionPoints.Add(point, path);
                        }
                    }
                }
            }


            //=========================================
            // Intrados
            //=========================================

            DataTree<Surface> intradosFaces = new DataTree<Surface>();

            foreach (GH_Path path in intersectionPoints.Paths)
            {
                var pts = intersectionPoints.Branch(path);
                if (pts != null && pts.Count == 4)
                {
                    // Ensure the points are in the correct order for surface creation
                    Point3d pt0 = pts[0];
                    Point3d pt1 = pts[1];
                    Point3d pt2 = pts[3];
                    Point3d pt3 = pts[2];
                    
                    // Create a 4-point surface
                    Surface srf = NurbsSurface.CreateFromCorners(pt0, pt1, pt2, pt3);
                    if (srf != null)
                        intradosFaces.Add(srf, path);
                        
                }
            }

            //=========================================
            // Average Normals
            //=========================================

            DataTree<Line> averageNormal = new DataTree<Line>();

            foreach (GH_Path path in intersectionCurve.Paths)
            {
                var line = intersectionCurve.Branch(path);
                if (line != null && line.Count == 4)
                {
                    //Average start
                    Point3d avStart = new Point3d(
                        line.Average(l => l.From.X),
                        line.Average(l => l.From.Y),
                        line.Average(l => l.From.Z)
                    );

                    //Average end
                    Point3d avEnd = new Point3d(
                        line.Average(l => l.To.X),
                        line.Average(l => l.To.Y),
                        line.Average(l => l.To.Z)
                    );

                    if (avStart.Z < avEnd.Z)
                    {
                        Line avLine = new Line(avStart, avEnd);
                        averageNormal.Add(avLine, path);
                    } 

                    if (avStart.Z > avEnd.Z)
                    {
                        Line avLine = new Line(avEnd, avStart);
                        averageNormal.Add(avLine, path);
                    }
                }
            }

            //=========================================
            // Extrados Planes
            //=========================================

            DataTree<Surface> extradosFaces = new DataTree<Surface>();

            foreach (GH_Path path in intradosPlanes.Paths)
            {
                // Get the original intrados plane
                var branch = intradosPlanes.get_Branch(path);
                if (branch == null || branch.Count == 0)
                    continue;

                Plane intradosPlane = ((GH_Plane)branch[0]).Value;

                // Get the corresponding average normal line
                var normalBranch = averageNormal.Branch(path);
                if (normalBranch == null || normalBranch.Count == 0)
                    continue;

                Line avgNormalLine = normalBranch[0];
                Vector3d direction = avgNormalLine.Direction;
                if (!direction.IsValid || direction.IsZero)
                    continue;

                direction.Unitize();

                // Ensure direction is opposite to intradosPlane.ZAxis
                if (direction * intradosPlane.ZAxis < 0)
                {
                    direction = -direction; // Reverse direction if not opposite
                }

                Vector3d moveVec = direction * thickness;

                // Move the plane origin
                Point3d newOrigin = intradosPlane.Origin + moveVec;
                Plane extradosPlane = new Plane(newOrigin, intradosPlane.XAxis, intradosPlane.YAxis);

                //=========================================
                // Intersect line with extrados plane
                //=========================================

                DataTree<Point3d> exIntersectionPoints = new DataTree<Point3d>();

                var linesBranch = intersectionCurve.Branch(path);
                if (linesBranch != null)
                {
                    foreach (var line in linesBranch)
                    {
                        if (Rhino.Geometry.Intersect.Intersection.LinePlane(line, extradosPlane, out double t))
                        {
                            Point3d point = line.PointAt(t);
                            exIntersectionPoints.Add(point, path);
                        }
                    }
                }
                //=========================================
                // Extrados Faces
                //=========================================

                foreach (GH_Path p in exIntersectionPoints.Paths)
                {
                    var pts = exIntersectionPoints.Branch(path);
                    if (pts != null && pts.Count == 4)
                    {
                        // Ensure the points are in the correct order for surface creation
                        Point3d pt0 = pts[2];
                        Point3d pt1 = pts[3];
                        Point3d pt2 = pts[1];
                        Point3d pt3 = pts[0];

                        // Create a 4-point surface
                        Surface srf = NurbsSurface.CreateFromCorners(pt0, pt1, pt2, pt3);
                        if (srf != null)
                            extradosFaces.Add(srf, path);

                    }
                }
            }

            //=========================================
            // Contact Faces
            //=========================================

            DataTree<Brep> contactFaces = new DataTree<Brep>();

            foreach (GH_Path path in intradosFaces.Paths)
            {
                Surface intrados = intradosFaces.Branch(path)[0];
                Brep intradosBrep = intrados.ToBrep();

                Surface extrados = extradosFaces.Branch(path)[0];
                Brep extradosBrep = extrados.ToBrep();

                // Extrair arestas
                List<Curve> edgesA = new List<Curve>();
                List<Curve> edgesB = new List<Curve>();

                for (int i = 0; i < intradosBrep.Edges.Count; i++)
                    edgesA.Add(intradosBrep.Edges[i].DuplicateCurve());

                for (int i = 0; i < extradosBrep.Edges.Count; i++)
                    edgesB.Add(extradosBrep.Edges[i].DuplicateCurve());

                // Reverse the list
                edgesB.Reverse();

                // Shift the list by one index forward (first element moves to the end)
                if (edgesB.Count > 1)
                {
                    Curve first = edgesB[0];
                    edgesB.RemoveAt(0);
                    edgesB.Add(first);
                }

                // Supondo que ambas têm o mesmo número de edges
                int edgeCount = Math.Min(edgesA.Count, edgesB.Count);

                for (int i = 0; i < edgeCount; i++)
                {
                    Curve c1 = edgesA[i];
                    Curve c2 = edgesB[i];
                    c2.Reverse();

                    Brep[] loft = Brep.CreateFromLoft(
                        new List<Curve> { c1, c2 },
                        Point3d.Unset,
                        Point3d.Unset,
                        LoftType.Straight,
                        false
                    );

                    if (loft != null && loft.Length > 0)
                    {
                        Brep loftedFace = loft[0];

                        // Adicionar ao DataTree
                        contactFaces.Add(loftedFace, path);
                    }
                }
            }

        //=========================================
        // Voussoirs
        //=========================================

        DataTree<Brep> voussoirs = new DataTree<Brep>();

            foreach (GH_Path path in contactFaces.Paths)
            {
                Brep contactFace0 = contactFaces.Branch(path)[0];
                Brep contactFace1 = contactFaces.Branch(path)[1];
                Brep contactFace2 = contactFaces.Branch(path)[2];
                Brep contactFace3 = contactFaces.Branch(path)[3];

                List<Brep> contactFaceList = new List<Brep> { contactFace0, contactFace1, contactFace2, contactFace3 };
                Brep[] jContactFaces = Brep.JoinBreps(contactFaceList, RhinoDoc.ActiveDoc.ModelAbsoluteTolerance);

                if (jContactFaces == null) continue;

                // Get the corresponding intrados and extrados faces
                Surface intrados = intradosFaces.Branch(path)[0];
                Surface extrados = extradosFaces.Branch(path)[0];

                Brep intradosBrep = intrados.ToBrep();
                Brep extradosBrep = extrados.ToBrep();

                List<Brep> brepsToJoin = new List<Brep> { intradosBrep, extradosBrep };
                brepsToJoin.AddRange(jContactFaces);

                Brep[] joined = Brep.JoinBreps(brepsToJoin, RhinoDoc.ActiveDoc.ModelAbsoluteTolerance);

                if (joined != null && joined.Length > 0)
                {
                    // Usually the first Brep is the solid
                    voussoirs.Add(joined[0], path);
                }
            }

            //log.Add("bye");

            DA.SetDataTree(0, voussoirs);
            DA.SetDataTree(2, intradosFaces);
            DA.SetDataTree(1, extradosFaces);
            DA.SetDataTree(3, contactFaces);
            //DA.SetDataTree(4, intersectionCurve);
        }
    }
}