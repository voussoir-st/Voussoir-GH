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
                              
                List<Brep> spBlocks = new List<Brep>();

                int spLineCounter = 0;
                foreach (Curve line in spLine)
                {
                    //System.Diagnostics.Debug.WriteLine($"Curve line in spLine");
                    // Create a vertical vector (Z direction) of length 10
                    Vector3d vertical = new Vector3d(0, 0, 10);

                    // Extrude the line along the vertical vector
                    Surface extrusion = Surface.CreateExtrusion(line, vertical);
                    Brep extrusionBrep = Brep.CreateFromSurface(extrusion);

                    GH_Structure<GH_Line> secondZline = new GH_Structure<GH_Line>();

                    //=========================================
                    // Voussoir interactions
                    //=========================================

                    //System.Diagnostics.Debug.WriteLine($""+ spLineCounter+"Just before the forloop. remappedVoussoirs.Branches.Count: " + remappedVoussoirs.Branches.Count);
                    //for (int branchIdx = 0; branchIdx < 3; branchIdx++)
                    for (int branchIdx = 0; branchIdx < remappedVoussoirs.Branches.Count; branchIdx++)
                    {
                            
                        System.Diagnostics.Debug.WriteLine($"remappedVoussoirs.Branches.Count " + remappedVoussoirs.Branches.Count);
                        var branch = remappedVoussoirs.Branches[branchIdx];
                        GH_Path branchPath = remappedVoussoirs.Paths[branchIdx];
                        GH_Structure<GH_Brep> intersectedVoussoirs = new GH_Structure<GH_Brep>();

                        int voussoirIndex = 0;
                        foreach (var vouss in branch)
                        {
                            Brep voussBrep = vouss.Value;

                            // Test for intersection (no boolean operation)
                            Curve[] intersectionCurves;
                            Point3d[] intersectionPoints;
                            bool intersects = Rhino.Geometry.Intersect.Intersection.BrepBrep(
                                voussBrep, extrusionBrep, RhinoDoc.ActiveDoc.ModelAbsoluteTolerance,
                                out intersectionCurves, out intersectionPoints
                            );

                            if (intersects && ((intersectionCurves != null && intersectionCurves.Length > 0) ||
                                               (intersectionPoints != null && intersectionPoints.Length > 0)))
                            {
                                //System.Diagnostics.Debug.WriteLine($"" + spLineCounter + " voussBrep: " + voussoirIndex);
                                //System.Diagnostics.Debug.WriteLine($""+spLineCounter+"INTERSECT got true. branchPath: "+ branchPath);

                                // Store the original voussoir in its respective branch
                                intersectedVoussoirs.Append(vouss, branchPath);

                                ////VIEW CONTENTS OF intersectedVoussoirs
                                //for (int i = 0; i < intersectedVoussoirs.Branches.Count; i++)
                                //{
                                //    GH_Path path = intersectedVoussoirs.Paths[i];
                                //    List<GH_Brep> branch2 = intersectedVoussoirs.Branches[i];

                                //    //System.Diagnostics.Debug.WriteLine($"Branch {i} at path {path}: {branch2.Count} items");

                                //    foreach (var ghBrep in branch2)
                                //    {
                                //        Brep brep = ghBrep?.Value;
                                //        if (brep == null) continue;

                                //        //System.Diagnostics.Debug.WriteLine($"  Brep: {brep.GetHashCode()} — Vertices: {brep.Vertices.Count}");
                                //    }
                                //}
                            }
                            voussoirIndex++;
                        }
                        int count2 = 0;
                        //System.Diagnostics.Debug.WriteLine($"" + spLineCounter + "intersectedVoussoirs.Branches.Count: " + intersectedVoussoirs.Branches.Count);

                        for (int i = 0; i < intersectedVoussoirs.Branches.Count; i++)
                        {
                            System.Diagnostics.Debug.WriteLine($"intersectedVoussoirs.Branches.Count " + intersectedVoussoirs.Branches.Count);
                            var b = intersectedVoussoirs.Branches[i];
                            GH_Path path = intersectedVoussoirs.Paths[i];

                            // Lista temporária para armazenar todos os valores Z deste ramo
                            List<double> zValues = new List<double>();
                            int a = 0;
                            foreach (var ghBrep in b)
                            {
                                //System.Diagnostics.Debug.WriteLine($"brep: {a}");

                                Brep brep = ghBrep.Value;
                                int q = 0;
                                // Iterar pelos vértices de cada Brep
                                foreach (BrepVertex v in brep.Vertices)
                                {
                                    
                                    zValues.Add(v.Location.Z);
                                    //System.Diagnostics.Debug.WriteLine($"vertex z: {q}");
                                    q++;
                                }
                                a++;
                            }

                            // Ordenar os valores de Z do maior para o menor
                            zValues.Sort();
                            zValues.Reverse();
                            Line zLine = new Line(spLine[0].PointAtStart, spLine[0].PointAtStart + Vector3d.ZAxis * zValues[1]);
                            secondZline.Append(new GH_Line(zLine), path);
                            //System.Diagnostics.Debug.WriteLine($"secondZline {count2}");
                            count2++;
                            //System.Diagnostics.Debug.WriteLine($"zValues:" + zValues);
                        }
                    }

                    Vector3d vert = line.PointAtStart + Vector3d.ZAxis * 1 - line.PointAtStart;
                    Line sLine = new Line(line.PointAtStart, line.PointAtEnd);
                    Rhino.Geometry.Plane spPlane = new Rhino.Geometry.Plane(line.PointAtStart, sLine.Direction, vert);
                    double d = spPlane.DistanceTo(sum);
                    //AddRuntimeMessage(GH_RuntimeMessageLevel.Remark, "debug 01 d: " + d);
                    
                    if (d > 0)
                    {
                        spPlane.Flip();
                        //AddRuntimeMessage(GH_RuntimeMessageLevel.Remark, "debug 02 d: " + d);
                    }
                    Vector3d offsetVec = spPlane.ZAxis * spWidth;
                    Line offsetLine = new Line(sLine.PointAt(0) + offsetVec, sLine.PointAt(1) + offsetVec);
                    offsetLine.Extend(20, 20);
                    List<Line> springerLines = new List<Line>();

                    foreach (GH_Plane ghPlane in tPlanes)
                    {
                        Rhino.Geometry.Plane p = ghPlane.Value;

                        Point3d ptA = Point3d.Unset;
                        Point3d ptB = Point3d.Unset;
                        bool hasA = false, hasB = false;

                        if (Rhino.Geometry.Intersect.Intersection.LinePlane(sLine, p, out double t))
                        {
                            ptA = sLine.PointAt(t);
                            hasA = true;
                        }

                        if (Rhino.Geometry.Intersect.Intersection.LinePlane(offsetLine, p, out double a))
                        {
                            ptB = offsetLine.PointAt(a);
                            hasB = true;
                        }

                        if (hasA && hasB)
                        {
                            Line aLine = new Line(ptA, ptB);
                            springerLines.Add(aLine);
                        }
                    }

                    //=========================================
                    // springer's base shape
                    //=========================================
                    Line bl = new Line(line.PointAtStart, line.PointAtEnd);

                    for (int i = 0; i < springerLines.Count - 1; i++)
                    {
                        Curve springerLine = springerLines[i].ToNurbsCurve();
                        Curve springerLinesNext = springerLines[i + 1].ToNurbsCurve();

                        Brep springerSurface = Brep.CreateFromLoft(new List<Curve> { springerLine, springerLinesNext }, Point3d.Unset, Point3d.Unset, LoftType.Normal, false)[0];
                        spBlocks.Add(springerSurface);
                        List<Curve> spSurfaceEdges = new List<Curve>();
                        Curve spSurfaceEdge0 = springerSurface.Edges[0];
                        spSurfaceEdges.Add(spSurfaceEdge0);
                        Curve spSurfaceEdge1 = springerSurface.Edges[1];
                        spSurfaceEdges.Add(spSurfaceEdge1);
                        Curve spSurfaceEdge2 = springerSurface.Edges[2];
                        spSurfaceEdges.Add(spSurfaceEdge2);
                        Curve spSurfaceEdge3 = springerSurface.Edges[3];
                        spSurfaceEdges.Add(spSurfaceEdge3);
                        Curve.JoinCurves(spSurfaceEdges, RhinoDoc.ActiveDoc.ModelAbsoluteTolerance);

                        Line l = secondZline.Branches[i][0].Value;
                        double a = l.Length;
                        Curve verticalEdge = new Line(
                            springerLine.PointAtStart,
                            springerLine.PointAtStart + Vector3d.ZAxis * a
                        ).ToNurbsCurve();

                        Brep[] springerBreps = Brep.CreateFromSweep(springerLine, spSurfaceEdges, true, RhinoDoc.ActiveDoc.ModelAbsoluteTolerance);
                        if (springerBreps != null && springerBreps.Length > 0)
                        {
                            Brep spb = springerBreps[0];
                            spBlocks.AddRange(springerBreps);
                        }
                        else
                        {
                            System.Diagnostics.Debug.WriteLine($"⚠️ Sweep failed at index {i}");
                        }
                    }
                    
                    spLineCounter++;
                    DA.SetDataTree(1, secondZline);
                }
                //DA.SetDataList(0, walls);
                //DA.SetDataTree(1, intersectedVoussoirs);
                //DA.SetDataTree(0, intersectedVoussoirs);
                //DA.SetDataTree(1, secondZline);
                DA.SetDataList(2, spBlocks);
            }
        }
    }
}