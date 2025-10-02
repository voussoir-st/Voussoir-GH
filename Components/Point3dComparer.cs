using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Rhino.Geometry;
using System.Collections.Generic;

namespace VoussoirPlugin03.Components
{
    public class Point3dComparer : IEqualityComparer<Point3d>
    {
        private double tolerance;

        public Point3dComparer(double tolerance = Rhino.RhinoMath.ZeroTolerance)
        {
            this.tolerance = tolerance;
        }

        public bool Equals(Point3d p1, Point3d p2)
        {
            return p1.DistanceTo(p2) <= tolerance;
        }

        public int GetHashCode(Point3d p)
        {
            unchecked
            {
                int hash = 17;
                hash = hash * 23 + p.X.GetHashCode();
                hash = hash * 23 + p.Y.GetHashCode();
                hash = hash * 23 + p.Z.GetHashCode();
                return hash;
            }
        }
    }
}