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
        public static class Utilsb
        {
            public static void OrientArcs(List<Curve> arcs)
            {
                if (arcs == null || arcs.Count != 2)
                    return;

                Point3d arc0Start = arcs[0].PointAtStart;
                Point3d arc0End = arcs[0].PointAtEnd;
                Point3d arc1Start = arcs[1].PointAtStart;
                Point3d arc1End = arcs[1].PointAtEnd;

                double dStart = arc0Start.DistanceTo(arc1Start);
                double dEnd = arc0End.DistanceTo(arc1End);
                double dCrossStartEnd = arc0Start.DistanceTo(arc1End);
                double dCrossEndStart = arc0End.DistanceTo(arc1Start);

                if (dStart + dEnd > dCrossStartEnd + dCrossEndStart)
                {
                    arcs[1].Reverse();
                }
            }
        }
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
                pManager.AddBrepParameter("Springers", "S", "Finished Springers", GH_ParamAccess.list);
                pManager.AddBrepParameter("Voussoirs", "V", "Non transformed voussoirs", GH_ParamAccess.tree);
                pManager.AddVectorParameter("Log", "L", "All messages generated during execution", GH_ParamAccess.list);
            }

            protected override void SolveInstance(IGH_DataAccess DA)
            {
                //=========================================
                // Declare and initialize input variables
                //=========================================
                List<Curve> spLine = new List<Curve>();
                GH_Structure<GH_Brep> voussoirs = new GH_Structure<GH_Brep>();
                List<Plane> tPlanes = new List<Plane>();
                double spWidth = 0.3;

                //=========================================
                // Retrieve input data from Grasshopper
                //=========================================
                if (!DA.GetDataList(0, spLine)) return;
                if (!DA.GetDataTree(1, out voussoirs)) return;
                if (!DA.GetDataList(2, tPlanes)) return;
                if (!DA.GetData(3, ref spWidth)) return;

                //=========================================
                // springer line intersection points
                //=========================================
                List<Point3d> splPtsA = new List<Point3d>();
                List<Point3d> splPtsB = new List<Point3d>();
                Utilsb.OrientArcs(spLine);
                List<Vector3d> sprVector = new List<Vector3d>();

                foreach (Plane pl in tPlanes)
                {
                    double t; // parameter along the line where intersection occurs
                    Line line = new Line(spLine[0].PointAtStart, spLine[0].PointAtEnd);
                    
                        if (Rhino.Geometry.Intersect.Intersection.LinePlane(line, pl, out t))
                        {
                            // Compute the 3D point at parameter t
                            Point3d pt = line.PointAt(t);
                            splPtsA.Add(pt);
                        }

                    double p; // parameter along the line where intersection occurs
                    Line lineB = new Line(spLine[1].PointAtStart, spLine[1].PointAtEnd);

                    if (Rhino.Geometry.Intersect.Intersection.LinePlane(lineB, pl, out p))
                    {
                        // Compute the 3D point at parameter t
                        Point3d pt = lineB.PointAt(p);
                        splPtsB.Add(pt);
                    }
                }
                for (int i = 0; i < splPtsA.Count; i++)
                {
                    Line spl = new Line(splPtsA[i], splPtsB[i]);
                    Vector3d splV = spl.PointAt(1) - spl.PointAt(0);
                    splV.Unitize();
                    sprVector.Add(splV);
                }
                DA.SetDataList(2, sprVector);
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

                //=========================================
                // Create Springers
                //=========================================
                List<Brep> springers = new List<Brep>();

                foreach (Curve line in spLine)
                {
                    // Create a vertical vector (Z direction) of length 10
                    Vector3d vertical = new Vector3d(0, 0, 10);

                    // Extrude the line along the vertical vector
                    Surface extrusion = Surface.CreateExtrusion(line, vertical);
                    Brep extrusionBrep = Brep.CreateFromSurface(extrusion);
                    
                    for (int branchIdx = 0; branchIdx < remappedVoussoirs.Branches.Count; branchIdx++)
                    {
                        var branch = remappedVoussoirs.Branches[branchIdx];
                        GH_Path branchPath = remappedVoussoirs.Paths[branchIdx];

                        GH_Structure<GH_Brep> intersectedVoussoirs = new GH_Structure<GH_Brep>();

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
                                // Store the original voussoir in its respective branch
                                intersectedVoussoirs.Append(vouss, branchPath);
                            }
                        }

                        GH_Structure<GH_Number> secondZValues = new GH_Structure<GH_Number>();

                        for (int i = 0; i < intersectedVoussoirs.Branches.Count; i++)
                        {
                            var b = intersectedVoussoirs.Branches[i];
                            GH_Path path = intersectedVoussoirs.Paths[i];

                            // Lista temporária para armazenar todos os valores Z deste ramo
                            List<double> zValues = new List<double>();

                            foreach (var ghBrep in b)
                            {
                                Brep brep = ghBrep.Value;

                                // Iterar pelos vértices de cada Brep
                                foreach (BrepVertex v in brep.Vertices)
                                {
                                    zValues.Add(v.Location.Z);
                                }
                            }

                            // Ordenar os valores de Z do maior para o menor
                            zValues.Sort();
                            zValues.Reverse();

                            // Se existirem pelo menos dois valores, guarda o de índice 1
                            if (zValues.Count > 1)
                            {
                                secondZValues.Append(new GH_Number(zValues[1]), path);
                            }
                            else if (zValues.Count == 1)
                            {
                                // Caso especial: só há um valor, guarda esse
                                secondZValues.Append(new GH_Number(zValues[0]), path);
                            }
                            else
                            {
                                // Caso o ramo esteja vazio, podes guardar null ou ignorar
                                secondZValues.Append(null, path);
                            }
                        }

                        //DA.SetDataTree(2, secondZValues);
                    }
                }
                //DA.SetDataList(0, walls);
                //DA.SetDataTree(1, intersectedVoussoirs);               
            }
        }
    }
}
