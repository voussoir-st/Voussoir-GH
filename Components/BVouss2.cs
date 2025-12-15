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
using System.Diagnostics;
using System.Linq;
using System.Security.Cryptography;
using VoussoirPlugin03.Properties;

namespace Components
{
    public class VoussoirCreate2 : GH_Component
    {
        public VoussoirCreate2()
         : base(
               "Create Voussoirs2",
               "BVouss2",
               "Creates Voussoirs based on predefined planes.",
               "Voussoir",
               "Core Geometry"
               )
        { }

        public override Guid ComponentGuid => new Guid("FC88D9F2-CD3B-4C41-ADFF-FD189794137C");

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
            pManager.AddTextParameter("Boundaries", "Boundaries", "Voussoir Boundaries (Indexes of Division Planes for each voussoir).", GH_ParamAccess.tree);
            pManager.AddNumberParameter("Voussoir Thickness", "Thickness", "User defined Voussoir Thickness", GH_ParamAccess.tree, 0.3);
        }

        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddBrepParameter("Voussoir", "V", "Voussoir Blocks", GH_ParamAccess.tree);
            pManager.AddBrepParameter("Extrados", "E", "Extrados Surfaces", GH_ParamAccess.tree);
            pManager.HideParameter(1);
            pManager.AddBrepParameter("Intrados", "I", "Intrados Surfaces", GH_ParamAccess.tree);
            pManager.HideParameter(2);
            pManager.AddBrepParameter("Contact Faces", "CF", "Voussoir contact faces", GH_ParamAccess.tree);
            pManager.HideParameter(3);
            //pManager.AddPlaneParameter("Log", "L", "All messages generated during execution", GH_ParamAccess.tree);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {

            //List<string> log = new List<string>();
            //log.Clear();
            GH_Structure<GH_Plane> intradosPlanesTree;
            GH_Structure<GH_Plane> divPlanesTree;
            GH_Structure<GH_String> divPlanesIndexesTree;
            GH_Structure<GH_Number> thicknesses;

            //log.Add("hello");

            //RhinoApp.WriteLine("SolveInstance called");

            if (!DA.GetDataTree(0, out intradosPlanesTree)) return;
            if (!DA.GetDataTree(1, out divPlanesTree)) return;
            if (!DA.GetDataTree(2, out divPlanesIndexesTree)) return;
            if (!DA.GetDataTree(3, out thicknesses)) return;

            double thickness = thicknesses[0][0].Value;

            //Per Vault
            GH_Structure<GH_Brep> voussoirsTree = new GH_Structure<GH_Brep>();
            GH_Structure<GH_Brep> extradosTree = new GH_Structure<GH_Brep>(); 
            GH_Structure<GH_Brep> intradosTree = new GH_Structure<GH_Brep>(); 
            GH_Structure<GH_Brep> contactFacesTree = new GH_Structure<GH_Brep>();
            GH_Structure<GH_Plane> log = new GH_Structure<GH_Plane>();

            foreach (GH_Path path in intradosPlanesTree.Paths)
            {
                if (thicknesses.Branches.Count > 1)
                    thickness = thicknesses[path][0].Value;
                List<Plane> intradosPlanes = intradosPlanesTree.Branches[intradosPlanesTree.Paths.IndexOf(path)].Select(x => x.Value).ToList();
                List<Plane> divPlanes = divPlanesTree.Branches[divPlanesTree.Paths.IndexOf(path)].Select(x => x.Value).ToList();
                List<String> divPlanesIndexes = divPlanesIndexesTree.Branches[divPlanesIndexesTree.Paths.IndexOf(path)].Select(x => x.Value).ToList();

                //Debug.WriteLine($"Path: {path}");

                //Per Intrados
                for (int p = 0; p < intradosPlanes.Count; p++)
                {
                    //Debug.WriteLine($"\nIntrados: {p}");
                    Plane intradosPlane = intradosPlanes[p];
                    String planeIndexString = divPlanesIndexes[p];
                    List<Plane> intradosDivPlanes = new List<Plane>();

                    //Get indexes for division planes
                    string trimmed = planeIndexString.TrimStart('B').Trim('{', '}');
                    string[] parts = trimmed.Split(';');
                    List<int> planeIndices = new List<int>();
                    foreach (string part in parts)
                    {
                        if (int.TryParse(part, out int index))
                        {
                            planeIndices.Add(index);
                        }
                    }
                    foreach (int index in planeIndices)
                    {
                        intradosDivPlanes.Add(divPlanes[index]);
                    }

                    //Debug.WriteLine($"intradosDivPlanes: {intradosDivPlanes.Count}");
                    //log.AppendRange(intradosDivPlanes.Select(pl => new GH_Plane(pl)), path);
                    //Plane Intersections                    
                    List<Line> intLines = new List<Line>();

                    //reorder planes
                    List<Point3d> origins = new List<Point3d>();
                    foreach (Plane pl in intradosDivPlanes)
                    {
                        origins.Add(pl.Origin);
                    }
                    //Debug.WriteLine($"origins: {origins.Count}");
                    Plane.FitPlaneToPoints(origins, out Plane fittedPlane);

                    var originsX = 0.0;
                    var originsY = 0.0;
                    var originsZ = 0.0;

                    foreach (Point3d q in origins)
                    {
                        originsX += q.X;
                        originsY += q.Y;
                        originsZ += q.Z;
                    }

                    var avgX = originsX / origins.Count;
                    var avgY = originsY / origins.Count;
                    var avgZ = originsZ / origins.Count;

                    var avgPoint = new Point3d(avgX, avgY, avgZ);

                    var fittedCircle = new Circle(fittedPlane, avgPoint, 1);

                    //planes.Append(new GH_Curve(fittedCircle.ToNurbsCurve()), path);
                    List<Plane> sortedPlanes = intradosDivPlanes.OrderBy(pl =>
                    {
                        double t;
                        fittedCircle.ClosestParameter(pl.Origin, out t);
                        return t;
                    }).ToList();

                    
                    //Debug.WriteLine($"sortedPlanes: {sortedPlanes.Count}");
                    //log.Append(new GH_Curve(fittedCircle.ToNurbsCurve()));

                    sortedPlanes.Add(sortedPlanes[0]);
                    for (int j = 0; j < sortedPlanes.Count - 1; j++)
                    {
                        Line line;
                        bool tf = Intersection.PlanePlane(sortedPlanes[j], sortedPlanes[j + 1], out line);
                        intLines.Add(line);
                    }
                    //planes.AppendRange(intLines.Select(pl => new GH_Line(pl)), path);

                    //new int lines
                    var newIntLines = new List<Line>();
                    foreach (Line l in intLines)
                    {
                        if (l.Direction * intradosPlane.ZAxis < 0)
                            l.Flip();
                        var t = thickness + 1;
                        var start = l.From;
                        var end = l.From + l.Direction * t;
                        newIntLines.Add(new Line(start, end));
                    }

                    //log.AppendRange(newIntLines.Select(l => new GH_Line(l)), path);

                    List<Point3d> intLinesStart = new List<Point3d>();
                    List<Point3d> intLinesEnd = new List<Point3d>();
                    foreach (Line l in newIntLines)
                    {
                        intLinesStart.Add(l.From);
                        intLinesEnd.Add(l.To);
                    }

                    //Debug.WriteLine($"intLines: {intLines.Count}");


                    //planes.AppendRange(newIntLines.Select(pl => new GH_Line(pl)), path);
                    //new normal
                    Point3d normalStart = Utils.AveragePoint(intLinesStart);
                    Point3d normalEnd = Utils.AveragePoint(intLinesEnd);
                    var newNormal = new Line(normalStart, normalEnd);
                    var newNormalVector = newNormal.Direction;
                    newNormalVector.Unitize();
                    if (newNormalVector * intradosPlane.ZAxis < 0)
                    {
                        newNormalVector = -newNormalVector; // Reverse direction if not opposite
                    }                                        

                    //Extrados Plane
                    var extradosPlaneOrigin = intradosPlane.Origin + newNormalVector * thickness;
                    var extradosPlane = new Plane(extradosPlaneOrigin, intradosPlane.XAxis, intradosPlane.YAxis);


                    //intrados Surface
                    List<Point3d> intradosPoints = new List<Point3d>();
                    //Debug.WriteLine($"intradosPoints: {intradosPoints.Count}");

                    int a = 0;
                    foreach (Line l in newIntLines)
                    {
                        //Debug.WriteLine($"l: {a}");
                        a++;
                        double t;
                        Intersection.LinePlane(l, intradosPlane, out t);
                        intradosPoints.Add(l.PointAt(t));
                    }
                    List<Point3d> intradosPolylinePoints = new List<Point3d>(intradosPoints);
                    intradosPolylinePoints.Add(intradosPoints[0]);
                    var intradosPolyline = new Polyline(intradosPolylinePoints);
                    //var intradosSurface = Brep.CreatePlanarBreps(new PolylineCurve(intradosPolyline), RhinoMath.ZeroTolerance);
                    var intradosSurface = NurbsSurface.CreateFromCorners(intradosPoints[0], intradosPoints[1], intradosPoints[2], intradosPoints[3]);
                    double u, v;
                    bool success = intradosSurface.ClosestPoint(intradosPlane.Origin, out u, out v);
                    var intradosNormal = intradosSurface.NormalAt(u, v);
                    intradosNormal.Unitize();
                    if (intradosNormal * intradosPlane.ZAxis > 0)
                        intradosSurface.Transpose();
                    intradosTree.Append(new GH_Brep(intradosSurface.ToBrep()), path);

                    //extrados surface
                    List<Point3d> extradosPoints = new List<Point3d>();
                    foreach (Line l in newIntLines)
                    {
                        double t;
                        Intersection.LinePlane(l, extradosPlane, out t);
                        extradosPoints.Add(l.PointAt(t));
                    }
                    List<Point3d> extradosPolylinePoints = new List<Point3d>(extradosPoints); 
                    extradosPolylinePoints.Add(extradosPoints[0]);
                    var extradosPolyline = new Polyline(extradosPolylinePoints);
                    //var extradosSurface = Brep.CreatePlanarBreps(new PolylineCurve(extradosPolyline), RhinoMath.ZeroTolerance);
                    var extradosSurface = NurbsSurface.CreateFromCorners(extradosPoints[0], extradosPoints[1], extradosPoints[2], extradosPoints[3]);
                    double o, n;
                    bool success1 = extradosSurface.ClosestPoint(intradosPlane.Origin, out o, out n);
                    var extradosNormal = extradosSurface.NormalAt(o, n);
                    extradosNormal.Unitize();
                    if (extradosNormal * intradosPlane.ZAxis < 0)
                        extradosSurface.Transpose();
                    extradosTree.Append(new GH_Brep(extradosSurface.ToBrep()), path);

                    //extradosTree.Append(new GH_Brep(extradosSurface[0].Faces[0].ToBrep()), path);
                    //Debug.WriteLine($"intradosPolylinePoints: {intradosPolylinePoints.Count}");
                    //Debug.WriteLine($"extradosPolylinePoints: {extradosPolylinePoints.Count}");
                    //Debug.WriteLine($"intradosPoints: {intradosPoints.Count}");
                    //Debug.WriteLine($"extradosPoints: {extradosPoints.Count}");

                    //contact faces
                    List<Brep> contactFaces = new List<Brep>();
                    for (int j = 0; j < intradosPoints.Count; j++)
                    {
                        List<Point3d> facePoints = new List<Point3d>();
                        facePoints.Add(intradosPolylinePoints[j]);
                        facePoints.Add(intradosPolylinePoints[j + 1]);
                        facePoints.Add(extradosPolylinePoints[j]);
                        facePoints.Add(extradosPolylinePoints[j + 1]);
                        Circle.TryFitCircleToPoints(facePoints, out Circle faceCircle);
                        List<Point3d> orderedpoints = facePoints.OrderBy(f =>
                        {
                            double t;
                            faceCircle.ClosestParameter(f, out t);
                            return t;
                        }).ToList();
                        var contactFaceSurface = NurbsSurface.CreateFromCorners(orderedpoints[0], orderedpoints[1], orderedpoints[2], orderedpoints[3]);
                        var contactFace = Brep.CreateFromSurface(contactFaceSurface);
                        contactFaces.Add(contactFace);
                    }
                    //Debug.WriteLine($"contactFaces: {contactFaces.Count}");
                    var joinedCF = Brep.JoinBreps(contactFaces, RhinoMath.ZeroTolerance);
                    contactFacesTree.Append(new GH_Brep(joinedCF[0]), path);
                    
                    //Voussoirs
                    List<Brep> voussoirParts = new List<Brep>();
                    voussoirParts.Add(intradosSurface.ToBrep());
                    voussoirParts.Add(extradosSurface.ToBrep());
                    voussoirParts.Add(joinedCF[0]);
                    var voussoir = Brep.JoinBreps(voussoirParts, RhinoDoc.ActiveDoc.ModelAbsoluteTolerance);                    
                    voussoirsTree.Append(new GH_Brep(voussoir[0]), path);                    
                }
            }

            DA.SetDataTree(0, voussoirsTree);
            DA.SetDataTree(1, extradosTree);
            DA.SetDataTree(2, intradosTree);
            DA.SetDataTree(3, contactFacesTree);
            //DA.SetDataTree(4, log);
        }
    }
}