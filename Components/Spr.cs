using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
using Rhino;
using Rhino.Geometry;
using Rhino.Render;
using System;
using System.Collections.Generic;


namespace VoussoirPlugin03
{
    namespace Components
    {
        public class VoussoirCreate : GH_Component
        {
            public VoussoirCreate()
             : base(
                   "Springer",
                   "Spr",
                   "Creates a Springer based on the voussoirs closest to the springer line",
                   "Voussoir",
                   "Vault Creation"
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
                pManager.AddCurveParameter("Springer Line", "SpringerLine", "List of base lines to create vault springers", GH_ParamAccess.list);
                pManager.AddBrepParameter("Voussoirs", "Voussoirs", "Voussoirs to analyse", GH_ParamAccess.tree);
                pManager.AddPlaneParameter("Transversal Planes", "TransversalPlanes", "Planes at each span division", GH_ParamAccess.list);
                pManager.AddNumberParameter("Springer Width", "SpringerWidth", "Distance perpendicular to springer line", GH_ParamAccess.item, 0.3);
            }

            protected override void RegisterOutputParams(GH_OutputParamManager pManager)
            {
                pManager.AddBrepParameter("Springers", "S", "Finished Springers", GH_ParamAccess.tree);
                pManager.AddLineParameter("Voussoirs", "V", "Non transformed voussoirs", GH_ParamAccess.tree);
                pManager.AddBrepParameter("Log", "L", "All messages generated during execution", GH_ParamAccess.list);
            }

            protected override void SolveInstance(IGH_DataAccess DA)
            {
                //=========================================
                // Declare and initialize input variables
                //=========================================
                List<Curve> spLine = new List<Curve>();
                GH_Structure<GH_Brep> voussoirs = new GH_Structure<GH_Brep>();
                List<GH_Plane> tPlanes = new List<GH_Plane>();
                double spWidth = 0.3;

                //=========================================
                // Retrieve input data from Grasshopper
                //=========================================
                if (!DA.GetDataList(0, spLine)) return;
                if (!DA.GetDataTree(1, out voussoirs)) return;
                if (!DA.GetDataList(2, tPlanes)) return;
                if (!DA.GetData(3, ref spWidth)) return;

                //=========================================
                // Remap voussoirs tree: {A;B} → {B}(A) (Path Mapper style)
                //=========================================
                GH_Structure<GH_Brep> remappedVoussoirs = new GH_Structure<GH_Brep>();

                for (int branchIdx = 0; branchIdx < voussoirs.PathCount; branchIdx++)
                {
                    GH_Path oldPath = voussoirs.Paths[branchIdx];
                    int[] oldIndices = oldPath.Indices;

                    // Only proceed if the path has at least two indices
                    if (oldIndices.Length >= 2)
                    {
                        int A = oldIndices[0];
                        int B = oldIndices[1];
                        var branch = voussoirs.Branches[branchIdx];

                        for (int itemIdx = 0; itemIdx < branch.Count; itemIdx++)
                        {
                            // New path is {B}
                            GH_Path newPath = new GH_Path(B);
                            // Item index in new branch is A
                            remappedVoussoirs.Insert(branch[itemIdx], newPath, A);
                        }
                    }
                    else
                    {
                        // If path does not have two indices, just append as is
                        var branch = voussoirs.Branches[branchIdx];
                        for (int itemIdx = 0; itemIdx < branch.Count; itemIdx++)
                        {
                            remappedVoussoirs.Append(branch[itemIdx], oldPath);
                        }
                    }
                }

                //=========================================
                // Extract extrados faces from voussoirs
                //=========================================
                GH_Structure<GH_Brep> extrados = new GH_Structure<GH_Brep>();

                foreach (GH_Brep b in remappedVoussoirs.AllData(true))
                {
                    Brep brep = b.Value;

                    if (brep.Faces.Count > 2)
                    {
                        BrepFace extradosFace = brep.Faces[1];
                        Brep extradosBrep = extradosFace.DuplicateFace(false);
                        extrados.Append(new GH_Brep(extradosBrep), new GH_Path(0));
                    }
                }
                
                int vertexCount = 0;
                double sumx = 0;
                double sumy = 0;
                double sumz = 0;
                
                foreach (GH_Brep ghBrep in remappedVoussoirs.AllData(true))
                {
                    Brep brep = ghBrep.Value;
                    if (brep == null) continue;
                    Point3d v = brep.Vertices[1].Location;                 
                    vertexCount += 1;
                    sumx += v.X;
                    sumy += v.Y;
                    sumz += v.Z;
                }
                double x = sumx / vertexCount;
                //System.Diagnostics.Debug.WriteLine($"x: {x}");
                double y = sumy / vertexCount;
                //System.Diagnostics.Debug.WriteLine($"y: {y}");
                double z = sumz / vertexCount;
                //System.Diagnostics.Debug.WriteLine($"z: {z}");
                Point3d sum = new Point3d(x, y, z);                

                //=========================================
                // Create Springers
                //=========================================
                List<Brep> springers = new List<Brep>();

                GH_Structure<GH_Brep> intersectedVoussoirs = new GH_Structure<GH_Brep>();
                GH_Structure<GH_Brep> springerVolumes = new GH_Structure<GH_Brep>();

                GH_Structure<GH_Line> secondZline = new GH_Structure<GH_Line>();
                var linecount = 0;
                //int spLineCounter = 0;

                foreach (Curve line in spLine)
                {
                    System.Diagnostics.Debug.WriteLine($"\n\nlinecount: " + linecount);
                    linecount++;
                    //=========================================
                    // Extrude line to create vertical reference
                    //=========================================
                    Vector3d vertical = Vector3d.ZAxis * 10;
                    Surface extrusion = Surface.CreateExtrusion(line, vertical);
                    Brep extrusionBrep = Brep.CreateFromSurface(extrusion);

                    // Estruturas locais por springer line
                    GH_Structure<GH_Line> secondZlineLocal = new GH_Structure<GH_Line>();
                    GH_Structure<GH_Brep> springerSurfacesLocal = new GH_Structure<GH_Brep>();
                    GH_Structure<GH_Brep> springerVolumesLocal = new GH_Structure<GH_Brep>();

                    //=========================================
                    // 1. Determinar interseções com voussoirs
                    //=========================================

                    for (int branchIdx = 0; branchIdx < remappedVoussoirs.Branches.Count; branchIdx++)
                    {
                        GH_Path branchPath = new GH_Path(branchIdx);
                        var branch = remappedVoussoirs.Branches[branchIdx];

                        foreach (var vouss in branch)
                        {
                            Brep voussBrep = vouss.Value;
                            if (voussBrep == null) continue;

                            if (Rhino.Geometry.Intersect.Intersection.BrepBrep(
                                voussBrep, extrusionBrep, RhinoDoc.ActiveDoc.ModelAbsoluteTolerance,
                                out Curve[] intersectionCurves, out Point3d[] intersectionPoints))
                            {
                                if ((intersectionCurves != null && intersectionCurves.Length > 0) ||
                                    (intersectionPoints != null && intersectionPoints.Length > 0))
                                {
                                    intersectedVoussoirs.Append(vouss, branchPath);
                                }
                            }
                        }
                    }

                    //=========================================
                    // 2. Calcular linhas Z (altura)
                    //=========================================
                    for (int i = 0; i < intersectedVoussoirs.Branches.Count; i++)
                    {
                        var branch = intersectedVoussoirs.Branches[i];
                        if (branch.Count == 0) continue;

                        List<double> zValues = new List<double>();
                        foreach (var ghBrep in branch)
                        {
                            Brep brep = ghBrep.Value;
                            foreach (BrepVertex v in brep.Vertices)
                                zValues.Add(v.Location.Z);
                        }

                        if (zValues.Count < 2) continue;
                        zValues.Sort(); zValues.Reverse();

                        double h = zValues[1];
                        Line zLine = new Line(line.PointAtStart, line.PointAtStart + Vector3d.ZAxis * h);
                        secondZlineLocal.Append(new GH_Line(zLine), new GH_Path(i));
                    }

                    //=========================================
                    // 3. Criar superfícies de base (springerSurfaces)
                    //=========================================
                    Vector3d vert = Vector3d.ZAxis;
                    Line sLine = new Line(line.PointAtStart, line.PointAtEnd);
                    Rhino.Geometry.Plane spPlane = new Rhino.Geometry.Plane(line.PointAtStart, sLine.Direction, vert);
                    double d = spPlane.DistanceTo(sum);
                    if (d > 0) spPlane.Flip();

                    Vector3d offsetVec = spPlane.ZAxis * spWidth;
                    Line offsetLine = new Line(sLine.PointAt(0) + offsetVec, sLine.PointAt(1) + offsetVec);
                    offsetLine.Extend(20, 20);

                    List<Line> springerLines = new List<Line>();
                    foreach (GH_Plane ghPlane in tPlanes)
                    {
                        Rhino.Geometry.Plane p = ghPlane.Value;
                        if (Rhino.Geometry.Intersect.Intersection.LinePlane(sLine, p, out double tA) &&
                            Rhino.Geometry.Intersect.Intersection.LinePlane(offsetLine, p, out double tB))
                        {
                            springerLines.Add(new Line(sLine.PointAt(tA), offsetLine.PointAt(tB)));
                        }
                    }

                    for (int i = 0; i < springerLines.Count - 1; i++)
                    {
                        Curve c0 = springerLines[i].ToNurbsCurve();
                        Curve c1 = springerLines[i + 1].ToNurbsCurve();
                        Brep[] lofts = Brep.CreateFromLoft(new List<Curve> { c0, c1 }, Point3d.Unset, Point3d.Unset, LoftType.Normal, false);
                        if (lofts != null && lofts.Length > 0)
                            springerSurfacesLocal.Append(new GH_Brep(lofts[0]), new GH_Path(i));
                    }
                    System.Diagnostics.Debug.WriteLine($"secondZlineLocal: " + secondZlineLocal.Branches.Count + "." + secondZlineLocal.Paths[0].Length + ";");
                    System.Diagnostics.Debug.WriteLine($"springerSurfacesLocal: " + springerSurfacesLocal.Paths.Count + ".");

                    //=========================================
                    // 4. Criar volumes a partir das superfícies e das linhas Z
                    //=========================================
                    var volumesran = 0;
                    List<Brep> spBlocks = new List<Brep>();

                    for (int s = 0; s < springerSurfacesLocal.Paths.Count; s++)
                    {
                        System.Diagnostics.Debug.WriteLine($"volumesran: " + volumesran + ".");

                        Brep surf = springerSurfacesLocal.Branches[s][0].Value;
                        Line hLine = secondZlineLocal.Branches[s][0].Value;
                        System.Diagnostics.Debug.WriteLine($"surf: " + surf.ObjectType + ";");
                        System.Diagnostics.Debug.WriteLine($"hLine: " + hLine.Length + ";");
                        double h = hLine.Length;
                        Vector3d moveVec = Vector3d.ZAxis * h;

                        Brep surf2 = surf.DuplicateBrep();
                        surf2.Translate(moveVec);
                        System.Diagnostics.Debug.WriteLine($"surf2: " + surf2.ObjectType + ";");

                        List<Curve> edges = new List<Curve>();
                        foreach (var edge in surf.Edges)
                        {
                            edges.Add(edge.EdgeCurve);
                        }

                        List<Curve> edges2 = new List<Curve>();
                        foreach (var edge in surf2.Edges)
                        {
                            edges2.Add(edge.EdgeCurve);
                        }

                        List<Brep> openBox = new List<Brep>();
                        for (int it = 0; it < edges.Count; it++)
                        {
                            var loftBreps = Brep.CreateFromLoft(new List<Curve> { edges[it], edges2[it] }, Point3d.Unset, Point3d.Unset, LoftType.Normal, false);
                            if (loftBreps != null && loftBreps.Length > 0)
                                openBox.Add(loftBreps[0]);
                        }
                        System.Diagnostics.Debug.WriteLine($"openBox: " + openBox.Count + ";");

                        List<Brep> box = new List<Brep>();
                        box.Add(surf);
                        box.Add(surf2);
                        box.AddRange(openBox);

                        System.Diagnostics.Debug.WriteLine($"box: " + box.Count + ";");
                        Brep[] joinedbox = Brep.JoinBreps(box, RhinoDoc.ActiveDoc.ModelAbsoluteTolerance);
                        spBlocks.Add(joinedbox[0]);
                        volumesran++;

                    }
                    System.Diagnostics.Debug.WriteLine($"spBlocks: " + spBlocks.Count + ";");
                }

                //DA.SetDataList(0, walls);
                DA.SetDataTree(0, intersectedVoussoirs);
                //DA.SetDataList(2, spBlocks);
                // DA.SetDataTree(1, secondZline);
            }
        }
    }
}