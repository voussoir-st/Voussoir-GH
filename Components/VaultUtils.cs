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
    }
}