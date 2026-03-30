using Grasshopper.Kernel;
using Grasshopper.Kernel.Parameters;
using Grasshopper.Kernel.Special;
using Rhino;
using Rhino.Geometry;
using Rhino.Geometry.Intersect;
using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;

namespace VoussoirPlugin03.Components.BaseSurface
{
    public class GroinVaultDefiningCurves : GH_Component
    {
        public GroinVaultDefiningCurves()
          : base(
            "Groin Vault Base Surface",
            "GVSrf",
            "Creates the base surface for the creation of a groin vault",
            "Voussoir",
            "Base Surface"
          )
        {
        }
        public override Guid ComponentGuid => new Guid("6335BBF4-5F06-4327-89F9-FFDE1C891A79");
        //protected override System.Drawing.Bitmap Icon
        //{
        //    get
        //    {
        //        //Use the imported Resources class directly
        //        return Resources.BVSrf;
        //    }
        //}
        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddCurveParameter("Springer Lines", "L1", "set(s) of 2 non intersecting lines", GH_ParamAccess.list);
            pManager.AddCurveParameter("Springer Lines", "L2", "set(s) of 2 non intersecting lines", GH_ParamAccess.list);
            pManager.AddNumberParameter("Vault Height", "H", "Arc's Height.", GH_ParamAccess.item, 2.0);
            //pManager.AddBooleanParameter("Span Direction", "SpanDirection", "0 = arcs on sides 0-1 and 2-3; 1 = arcs on sides 1-2 and 3-0.", GH_ParamAccess.item, true);
            pManager.AddIntegerParameter("Arc Profile", "A",
                "Right-click to select type of curve\n\nArc = 0\nParabola = 1\nCatenary = 2\n\n... or input one of the above integers", GH_ParamAccess.item, 0);

            var param = (Param_Integer)pManager[3];
            param.AddNamedValue("Arc", 0);
            param.AddNamedValue("Parabola", 1);
            param.AddNamedValue("Catenary", 2);

            pManager[1].Optional = true; // Vault Height
            pManager[2].Optional = true; // Span Direction
        }
        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddBrepParameter("Vault Surface", "S", "The lofted base surface between the two arcs.", GH_ParamAccess.tree);

            pManager.AddCurveParameter("Springer Lines", "L", "The 2 Horizontal lines (remaining sides).", GH_ParamAccess.tree);
            pManager.HideParameter(1);
            pManager.AddCurveParameter("Vault Arcs", "A", "The 2 generated arcs.", GH_ParamAccess.tree);
            pManager.HideParameter(2);
        }
        private List<Point3d> _previewPts = new List<Point3d>();

        public override bool IsPreviewCapable => true;

        public override void DrawViewportWires(IGH_PreviewArgs args)
        {
            base.DrawViewportWires(args);

            foreach (var pt in _previewPts)
            {
                args.Display.DrawPoint(pt, Rhino.Display.PointStyle.RoundSimple, 3, System.Drawing.Color.DarkRed);
            }
        }
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            double tol = 0.0001;
            //_previewPts.Clear();
            // Declare input variables
            List<Curve> springerLines = new List<Curve>();
            List<Curve> springerLines2 = new List<Curve>();
            double height = 2.0;
            bool spanDir = false;
            int mode = 0;

            // Ler inputs
            DA.GetDataList(0, springerLines);
            DA.GetDataList(1, springerLines2);
            DA.GetData(2, ref height);
            DA.GetData(3, ref mode);

            Curve c0 = springerLines[0];
            Curve c1 = springerLines[1];
            Curve c2 = springerLines2[0];
            Curve c3 = springerLines2[1];

            var int02 = CurvesIntersect(c0, c2);
            var int03 = CurvesIntersect(c0, c3);
            var int12 = CurvesIntersect(c1, c2);
            var int13 = CurvesIntersect(c1, c3);

            bool coplanar = int02.boolean && int03.boolean && int12.boolean && int13.boolean;

            if (!coplanar)
            {
                if (!coplanar)
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "The four Springer Lines must be coplanar.");
                    return;
                }
            }

            // Processar os dois conjuntos
            var result1 = ProcessSpringerPair(springerLines, height, spanDir, mode);
            var result2 = ProcessSpringerPair(springerLines2, height, spanDir, mode);

            List<Surface> surfaces = new List<Surface>();
            surfaces.Add(result1.loft);
            surfaces.Add(result2.loft);

            List<Curve> lines = new List<Curve>();
            lines.AddRange(result1.lines);
            lines.AddRange(result2.lines);

            List<Curve> arcs = new List<Curve>();
            arcs.AddRange(result1.arcs);
            arcs.AddRange(result2.arcs);

            List<Point3d> intPoints = new List<Point3d>();
            intPoints.Add(int02.points[0].PointA);
            intPoints.Add(int03.points[0].PointA);
            intPoints.Add(int12.points[0].PointA);
            intPoints.Add(int13.points[0].PointA);

            Plane basePlane;
            double er;
            var ok = Plane.FitPlaneToPoints(intPoints, out basePlane, out er);

            if (basePlane.ZAxis * Plane.WorldXY.ZAxis < 0) basePlane.Flip();

            var arcMidPoints = arcs
                .Where(c => c != null && c.IsValid)
                .Select(c => c.PointAtNormalizedLength(0.5))
                .OrderBy(p => basePlane.DistanceTo(p))     // ordena do menor Z para o maior
                .ToList();

            Point3d lowestPoint = arcMidPoints.First();   // só o mais baixo

            var newHeight = basePlane.DistanceTo(lowestPoint);

            var result3 = ProcessSpringerPair(springerLines, newHeight, spanDir, mode);
            var result4 = ProcessSpringerPair(springerLines2, newHeight, spanDir, mode);

            List<Surface> surfaces2 = new List<Surface>();
            surfaces2.Add(result3.loft);
            surfaces2.Add(result4.loft);

            List<Curve> lines2 = new List<Curve>();
            if (result3.lines[0].PointAtStart.DistanceTo(result3.lines[1].PointAtEnd) < result3.lines[0].PointAtStart.DistanceTo(result3.lines[1].PointAtStart)) ;
            {
                result3.lines[1].Reverse();
            }
            lines2.AddRange(result3.lines);
            if (result4.lines[0].PointAtStart.DistanceTo(result4.lines[1].PointAtEnd) < result4.lines[0].PointAtStart.DistanceTo(result4.lines[1].PointAtStart)) ;
            {
                result4.lines[1].Reverse();
            }
            lines2.AddRange(result4.lines);

            List<Curve> arcs2 = new List<Curve>();
            arcs2.AddRange(result3.arcs);
            arcs2.AddRange(result4.arcs);

            var arcMidPoints2 = arcs2
                .Where(c => c != null && c.IsValid)
                .Select(c => c.PointAtNormalizedLength(0.5))
                .OrderBy(p => basePlane.DistanceTo(p))     // ordena do menor Z para o maior
                .ToList();

            var highestZ2 = basePlane.DistanceTo(arcMidPoints2.Last());
            var lowestZ2 = basePlane.DistanceTo(arcMidPoints2.First());
            //System.Diagnostics.Debug.WriteLine($"---- highestZ2: {highestZ2}");
            //System.Diagnostics.Debug.WriteLine($"---- lowestZ2: {lowestZ2}");

            var shrinkFactor = lowestZ2 / highestZ2;

            //System.Diagnostics.Debug.WriteLine($"---- shrink: {shrinkFactor}");
            double srf0_height = Math.Max(
                basePlane.DistanceTo(arcMidPoints2[0]),
                basePlane.DistanceTo(arcMidPoints2[1])
            );

            double srf1_height = Math.Max(
                basePlane.DistanceTo(arcMidPoints2[2]),
                basePlane.DistanceTo(arcMidPoints2[3])
            );

            // Determine highest and lowest surfaces
            Surface highestSurface;
            Surface lowestSurface;

            if (srf0_height > srf1_height)
            {
                highestSurface = surfaces2[0];
                lowestSurface = surfaces2[1];
            }
            else
            {
                highestSurface = surfaces2[1];
                lowestSurface = surfaces2[0];
            }

            Surface scaledHigh = lowestSurface.Duplicate() as Surface;

            Transform scale1D = Transform.Scale(basePlane, 1, 1, shrinkFactor);

            scaledHigh.Transform(scale1D);

            var SAarcLine1 = new Line(result3.lines[0].PointAtStart, result3.lines[1].PointAtStart);
            var SA1midpoint = SAarcLine1.PointAt(0.5);

            var SAarcLine2 = new Line(result3.lines[0].PointAtEnd, result3.lines[1].PointAtEnd);
            var SA2midpoint = SAarcLine2.PointAt(0.5);

            var SBarcLine1 = new Line(result4.lines[0].PointAtStart, result4.lines[1].PointAtStart);
            var SB1midpoint = SBarcLine1.PointAt(0.5);

            var SBarcLine2 = new Line(result4.lines[0].PointAtEnd, result4.lines[1].PointAtEnd);
            var SB2midpoint = SBarcLine2.PointAt(0.5);

            Curve[] intCurves;
            Point3d[] interPoints;           
            Intersection.SurfaceSurface(scaledHigh, highestSurface, tol, out intCurves, out interPoints);
            var ic = intCurves.ToList();
            //System.Diagnostics.Debug.WriteLine("---- DEBUG: Surface Split ----");

            //System.Diagnostics.Debug.WriteLine("Intersection curves: " + intCurves.Length);

            Brep brepA = scaledHigh.ToBrep();
            Brep brepB = highestSurface.ToBrep();

            List<Curve> intCurvesList = new List<Curve>();

            List<Brep> breps = new List<Brep>();
            breps.Add(brepA);
            breps.Add(brepB);

            //_previewPts.Add(new Point3d(0, 0, 0));

            if (intCurves.Length > 2)
            {
                //System.Diagnostics.Debug.WriteLine("if");
                var curve0 = intCurves[0];
                var curve1 = intCurves[1];
                var curve2 = intCurves[2];
                var curve3 = intCurves[3];

                var highestPoint =
                    new[] { curve0, curve1, curve2, curve3 }
                    .SelectMany(c => new[] { c.PointAtStart, c.PointAtEnd })
                    .OrderByDescending(p => p.Z)
                    .First();

                curve0 = EnsureEndClosestToApex(curve0, highestPoint);
                curve1 = EnsureEndClosestToApex(curve1, highestPoint);
                curve2 = EnsureEndClosestToApex(curve2, highestPoint);
                curve3 = EnsureEndClosestToApex(curve3, highestPoint);

                List<Point3d> startpoints = new List<Point3d>();
                startpoints.Add(curve0.PointAtStart);
                startpoints.Add(curve1.PointAtStart);
                startpoints.Add(curve2.PointAtStart);
                startpoints.Add(curve3.PointAtStart);

                var startavg = new Point3d(
                    startpoints.Average(p => p.X),
                    startpoints.Average(p => p.Y),
                    startpoints.Average(p => p.Z)
                );

                var avgRadius = startavg.DistanceTo(curve0.PointAtStart);

                var orderingCircle = new Circle(startavg, avgRadius);

                List<Curve> newCurves = new List<Curve>();
                newCurves.Add(curve0);
                newCurves.Add(curve1);
                newCurves.Add(curve2);
                newCurves.Add(curve3);

                intCurves = newCurves
                    .Select(crv =>
                    {
                        // representative point on the curve
                        Point3d pt = crv.PointAtStart;

                        // find closest point on the ordering circle
                        double t;
                        orderingCircle.ClosestParameter(pt, out t);
                        //System.Diagnostics.Debug.WriteLine($"Curve start: {pt}, T = {t}");

                        return new { Curve = crv, T = t };
                    })
                    .OrderBy(x => x.T)
                    .Select(x => x.Curve)
                    .ToArray();

                Curve[] join01 = Curve.JoinCurves(new Curve[] { intCurves[0], intCurves[2] }, tol);

                Curve joined01 = join01.Length > 0 ? join01[0] : null;

                Curve[] join02 = Curve.JoinCurves(new Curve[] { intCurves[1], intCurves[3] }, tol);

                Curve joined02 = join02.Length > 0 ? join02[0] : null;


                intCurvesList.Add(joined01);
                intCurvesList.Add(joined02);
            }
            else
            {
                //System.Diagnostics.Debug.WriteLine("else");

                Curve curveA = intCurves[0];
                Curve curveB = intCurves[1];

                // 1. Encontrar o ponto mais próximo entre as duas curvas (o ápice do V)
                Point3d pa, pb;
                curveA.ClosestPoints(curveB, out pa, out pb);
                var dist = pa.DistanceTo(pb);
                // 2. Converter esses pontos em parâmetros t
                double ta, tb;
                curveA.ClosestPoint(pa, out ta);
                curveB.ClosestPoint(pb, out tb);

                // 3. Garantir que os parâmetros estão dentro do domínio
                Interval domA = curveA.Domain;
                Interval domB = curveB.Domain;

                ta = Math.Max(domA.Min, Math.Min(domA.Max, ta));
                tb = Math.Max(domB.Min, Math.Min(domB.Max, tb));

                // 4. Dividir as curvas no ápice
                Curve[] splitA = curveA.Split(ta);
                Curve[] splitB = curveB.Split(tb);

                var curve0 = splitA[0];
                var curve1 = splitA[1];
                var curve2 = splitB[0];
                var curve3 = splitB[1];

                curve0 = ShortenEnd(curve0, dist * 2);
                curve1 = ShortenStart(curve1, dist * 2);
                curve2 = ShortenEnd(curve2, dist * 2);
                curve3 = ShortenStart(curve3, dist * 2);

                curve0 = EnsureEndClosestToApex(curve0, pa);
                curve1 = EnsureEndClosestToApex(curve1, pa);
                curve2 = EnsureEndClosestToApex(curve2, pa);
                curve3 = EnsureEndClosestToApex(curve3, pa);

                List<Point3d> startpoints = new List<Point3d>();
                startpoints.Add(curve0.PointAtStart);
                startpoints.Add(curve1.PointAtStart);
                startpoints.Add(curve2.PointAtStart);
                startpoints.Add(curve3.PointAtStart);

                var startavg = new Point3d(
                    startpoints.Average(p => p.X),
                    startpoints.Average(p => p.Y),
                    startpoints.Average(p => p.Z)
                );

                var avgRadius = startavg.DistanceTo(curve0.PointAtStart);

                var orderingCircle = new Circle(startavg, avgRadius);

                List<Curve> newCurves = new List<Curve>();
                newCurves.Add(curve0);
                newCurves.Add(curve1);
                newCurves.Add(curve2);
                newCurves.Add(curve3);

                intCurves = newCurves
                    .Select(crv =>
                    {
                        // representative point on the curve
                        Point3d pt = crv.PointAtStart;

                        // find closest point on the ordering circle
                        double t;
                        orderingCircle.ClosestParameter(pt, out t);
                        //System.Diagnostics.Debug.WriteLine($"Curve start: {pt}, T = {t}");

                        return new { Curve = crv, T = t };
                    })
                    .OrderBy(x => x.T)
                    .Select(x => x.Curve)
                    .ToArray();

                List<Point3d> endpoints = new List<Point3d>();
                endpoints.Add(curve0.PointAtEnd);
                endpoints.Add(curve1.PointAtEnd);
                endpoints.Add(curve2.PointAtEnd);
                endpoints.Add(curve3.PointAtEnd);

                var endavg = new Point3d(
                    endpoints.Average(p => p.X),
                    endpoints.Average(p => p.Y),
                    endpoints.Average(p => p.Z)
                );

                var movedAvg = new Point3d(endavg.X, endavg.Y, endavg.Z + newHeight);
                var pullSrfAvg = brepA.ClosestPoint(movedAvg);

                Point3d[] intCurve1Points;
                var divIntCurve1 = intCurves[0].DivideByCount(10, true, out intCurve1Points);

                Point3d[] intCurve2Points;
                var divIntCurve2 = intCurves[1].DivideByCount(10, true, out intCurve2Points);

                intCurves[2].Reverse();
                Point3d[] intCurve3Points;
                var divIntCurve3 = intCurves[2].DivideByCount(10, true, out intCurve3Points);
                intCurves[3].Reverse();
                Point3d[] intCurve4Points;
                var divIntCurve4 = intCurves[3].DivideByCount(10, true, out intCurve4Points);

                List<Point3d> divCrvAPts = new List<Point3d>();
                divCrvAPts.AddRange(intCurve1Points);
                divCrvAPts.Add(pullSrfAvg);
                divCrvAPts.AddRange(intCurve3Points);

                //_previewPts.AddRange(divCrvAPts);

                var divCrvA = Curve.CreateInterpolatedCurve(divCrvAPts, 3);

                List<Point3d> divCrvBPts = new List<Point3d>();
                divCrvBPts.AddRange(intCurve2Points);
                divCrvBPts.Add(pullSrfAvg);
                divCrvBPts.AddRange(intCurve4Points);

                var divCrvB = Curve.CreateInterpolatedCurve(divCrvBPts, 3);

                intCurvesList.Add(divCrvA);
                intCurvesList.Add(divCrvB);
            }
            
            List<Curve> pullcrvs1 = new List<Curve>();
            foreach (Curve c in intCurvesList)
            {
                c.Extend(height, height);
                c.PullToBrepFace(brepA.Faces[0], tol);
                pullcrvs1.Add(c);
            }
            List<Curve> pullcrvs2 = new List<Curve>();
            foreach (Curve c in intCurvesList)
            {
                c.Extend(height, height);
                c.PullToBrepFace(brepB.Faces[0], tol);
                pullcrvs2.Add(c);
            }

            var splitSurfaceA = brepA.Split(pullcrvs1, tol);

            var splitSurfaceB = brepB.Split(pullcrvs2, tol);

            var SA1closestBrep = splitSurfaceA
                .OrderBy(b => b.ClosestPoint(SA1midpoint).DistanceTo(SA1midpoint))
                .First();

            var SA2closestBrep = splitSurfaceA
                .OrderBy(b => b.ClosestPoint(SA2midpoint).DistanceTo(SA2midpoint))
                .First();

            var SB1closestBrep = splitSurfaceB
                .OrderBy(b => b.ClosestPoint(SB1midpoint).DistanceTo(SB1midpoint))
                .First();

            var SB2closestBrep = splitSurfaceB
                .OrderBy(b => b.ClosestPoint(SB2midpoint).DistanceTo(SB2midpoint))
                .First();

            List<Brep> trimmedSurfaces = new List<Brep>();
            trimmedSurfaces.Add(SA1closestBrep);
            trimmedSurfaces.Add(SA2closestBrep);
            trimmedSurfaces.Add(SB1closestBrep);
            trimmedSurfaces.Add(SB2closestBrep);

            // Output das superfícies
            DA.SetDataList(0, trimmedSurfaces);
            DA.SetDataList(1, lines2);
            DA.SetDataList(2, ic);
        }
        Curve EnsureEndClosestToApex(Curve c, Point3d apex)
        {
            double dStart = c.PointAtStart.DistanceTo(apex);
            double dEnd = c.PointAtEnd.DistanceTo(apex);

            // If start is closer than end → flip
            if (dStart < dEnd)
                c.Reverse();

            return c;
        }
        Curve ShortenStart(Curve c, double length)
        {
            double t0 = c.Domain.Min;
            double t1 = c.Domain.Max;

            double newStart = t0 + length;
            if (newStart >= t1) return c;

            return c.Trim(new Interval(newStart, t1));
        }
        Curve ShortenEnd(Curve c, double length)
        {
            double t0 = c.Domain.Min;
            double t1 = c.Domain.Max;

            double newEnd = t1 - length;
            if (newEnd <= t0) return c; // avoid invalid trim

            return c.Trim(new Interval(t0, newEnd));
        }
        private (bool boolean, CurveIntersections points) CurvesIntersect(Curve a, Curve b)
        {
            var points = Intersection.CurveCurve(a, b, 0.001, 0.001);
            bool boolean = points != null && points.Count > 0;
            return (boolean, points);
        }
        private (bool ok, Surface loft, List<Curve> lines, List<Curve> arcs)
        ProcessSpringerPair(List<Curve> springerLines, double height, bool spanDir, int mode)
        {
            // Validar
            if (springerLines == null || springerLines.Count < 2)
                return (false, null, null, null);

            Utils.Utils.OrientArcs(springerLines);

            // Verificar interseção
            var events = Rhino.Geometry.Intersect.Intersection.CurveCurve(
                springerLines[0], springerLines[1],
                RhinoDoc.ActiveDoc.ModelAbsoluteTolerance, 0.0);

            if (events != null && events.Count > 0)
                return (false, null, null, null);

            // Definir edges
            var edges = new (Point3d A, Point3d B)[]
            {
                (springerLines[0].PointAtStart, springerLines[1].PointAtStart),
                (springerLines[1].PointAtStart, springerLines[1].PointAtEnd),
                (springerLines[1].PointAtEnd,   springerLines[0].PointAtEnd),
                (springerLines[0].PointAtEnd,   springerLines[0].PointAtStart)
            };

            // Escolher lados dos arcos
            int[][] arcPairs = new int[][]
            {
                new[] { 0, 2 },
                new[] { 1, 3 }
            };

            int sd = spanDir ? 1 : 0;
            var arcIdx = arcPairs[sd];
            var lineIdx = sd == 0 ? new[] { 1, 3 } : new[] { 0, 2 };

            // Clamp da altura
            if (mode == 0)
            {
                double dist0 = edges[arcIdx[0]].A.DistanceTo(edges[arcIdx[0]].B);
                double dist1 = edges[arcIdx[1]].A.DistanceTo(edges[arcIdx[1]].B);
                double minDist = Math.Min(dist0, dist1);
                if (height > minDist / 2)
                    height = minDist / 2;
            }

            // Criar arcos
            var arcs = new List<Curve>();
            foreach (int ei in arcIdx)
            {
                var (A, B) = edges[ei];
                var c = BuildSpanCurve(A, B, height, mode);
                //Curve d = c.Rebuild(150, 3, true);
                if (c != null) arcs.Add(c);
            }

            Utils.Utils.OrientArcs(arcs);

            // Criar linhas
            var lines = new List<Curve>();
            foreach (int ei in lineIdx)
            {
                var (A, B) = edges[ei];
                lines.Add(new LineCurve(A, B));
            }

            // Loft
            var loftBreps = Brep.CreateFromLoft(arcs, Point3d.Unset, Point3d.Unset, LoftType.Normal, false);
            if (loftBreps == null || loftBreps.Length == 0)
                return (false, null, lines, arcs);

            var alignedSurface = Components.BaseSurface.PolylineUtils.AlignNormalToWorldZ(loftBreps[0].Surfaces[0]);

            return (true, alignedSurface, lines, arcs);
        }
        private Curve BuildSpanCurve(Point3d start, Point3d end, double h, int mode)
        {
            if (h <= RhinoMath.ZeroTolerance)
                return new LineCurve(start, end); // altura nula → linha

            switch (mode)
            {
                case 0:
                    // Arco 3 pontos
                    {
                        double dist = start.DistanceTo(end);
                        if (h > dist)
                            h = dist; // Clamp height to the arc's span

                        Point3d mid = 0.5 * (start + end) + Vector3d.ZAxis * h;
                        var arc = new Arc(start, mid, end);
                        if (arc.IsValid) return arc.ToNurbsCurve();
                        return NurbsCurve.CreateInterpolatedCurve(new[] { start, mid, end }, 3);
                    }
                case 2:
                    // Catenária com extremos a 0 e meio a +h
                    return BuildCatenary(start, end, h, 100);
                case 1:
                default:
                    // Parábola por 3 pontos
                    {
                        Point3d mid = 0.5 * (start + end) + Vector3d.ZAxis * h;
                        return NurbsCurve.CreateInterpolatedCurve(new[] { start, mid, end }, 3);
                    }
            }
        }
        private Curve BuildCatenary(Point3d start, Point3d end, double rise, int samples)
        {
            // Draw catenary in a vertical plane between start and end (Z axis is up)
            if (samples < 3) samples = 3;

            // Project start and end to XY, but keep their Z values for later
            var startXY = new Point2d(start.X, start.Y);
            var endXY = new Point2d(end.X, end.Y);
            Vector2d dirXY = endXY - startXY;
            double L = dirXY.Length;
            if (L < RhinoMath.ZeroTolerance)
                return new LineCurve(start, end);
            dirXY.Unitize();

            // Find midpoint in XY
            var midXY = 0.5 * (startXY + endXY);

            // Solve catenary parameter 'a'
            double a = SolveCatenaryA(L, rise);
            if (double.IsNaN(a) || a <= 0)
            {
                // fallback: parabola in vertical plane
                Point3d apex = new Point3d(midXY.X, midXY.Y, Math.Max(start.Z, end.Z) + rise);
                return NurbsCurve.CreateInterpolatedCurve(new[] { start, apex, end }, 2);
            }

            var pts = new List<Point3d>(samples);
            for (int i = 0; i < samples; i++)
            {
                double t = i / (double)(samples - 1); // 0..1
                double x = (t - 0.5) * L;             // -L/2 .. L/2
                double z = rise - (a * Math.Cosh(x / a) - a); // vertical displacement
                // Interpolate XY position
                double px = midXY.X + x * dirXY.X;
                double py = midXY.Y + x * dirXY.Y;
                // Interpolate Z between start and end, then add catenary rise
                double baseZ = start.Z + t * (end.Z - start.Z);
                pts.Add(new Point3d(px, py, baseZ + z));
            }
            // Set exact start/end points
            if (pts.Count >= 2)
            {
                pts[0] = start;
                pts[pts.Count - 1] = end;
            }
            return Curve.CreateInterpolatedCurve(pts, 3);
        }
        private double SolveCatenaryA(double span, double rise)
        {
            // Resolve f(a) = a*cosh((span/2)/a) - a - rise = 0
            if (rise <= 0 || span <= RhinoMath.ZeroTolerance) return double.NaN;
            double a = Math.Max(span / 4.0, 1e-4); // palpite inicial conservador
            const double tol = 1e-10;
            const int maxIt = 100;
            for (int i = 0; i < maxIt; i++)
            {
                double x = (span * 0.5) / a;
                double cosh = Math.Cosh(x);
                double sinh = Math.Sinh(x);
                double f = a * cosh - a - rise;
                double df = cosh - 1.0 - (span * 0.5) * sinh / a; // derivada
                if (Math.Abs(df) < 1e-14) { a *= 1.2; continue; }
                double aNext = a - f / df;
                if (Math.Abs(aNext - a) < tol) return aNext;
                a = Math.Max(aNext, 1e-6);
            }
            return a; // devolve melhor estimativa
        }
        public override void AddedToDocument(GH_Document document)
        {
            base.AddedToDocument(document);

            if (Params.Input.Count < 3)
                return;

            Param_Integer param = Params.Input[3] as Param_Integer;
            if (param == null || param.SourceCount > 0)
                return;

            Attributes.PerformLayout();

            int x = (int)param.Attributes.Pivot.X - 135;
            int y = (int)param.Attributes.Pivot.Y - 11;

            GH_ValueList valueList = new GH_ValueList();
            valueList.CreateAttributes();

            valueList.ListItems.Clear();
            valueList.ListItems.Add(new GH_ValueListItem("Arc", "0"));
            valueList.ListItems.Add(new GH_ValueListItem("Parabola", "1"));
            valueList.ListItems.Add(new GH_ValueListItem("Catenary", "2"));

            valueList.SelectItem(0);

            valueList.Attributes.Pivot = new PointF(x, y);
            valueList.Attributes.ExpireLayout();

            document.AddObject(valueList, false);

            param.AddSource(valueList);
        }
    }
}
