using Rhino.Geometry;
using System.Collections.Generic;

namespace Components
{
    public static class Utils
    {
        /// <summary>
        /// Ensures that two curves (arcs) in the input list are oriented in the same parameter direction.
        /// If not, reverses the second arc in-place.
        /// </summary>
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
        public static Point3d AveragePoint(List<Point3d> pts)
        {
            if (pts == null || pts.Count == 0)
                return Point3d.Unset;

            double x = 0, y = 0, z = 0;

            foreach (var p in pts)
            {
                x += p.X;
                y += p.Y;
                z += p.Z;
            }

            double n = pts.Count;
            return new Point3d(x / n, y / n, z / n);
        }
    }
}