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
                //pManager.AddCurveParameter("Log", "L", "All messages generated during execution", GH_ParamAccess.tree);
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
                    }
                }
                //DA.SetDataList(0, walls);
                //DA.SetDataTree(1, intersectedVoussoirs);               
            }
        }
    }
}
