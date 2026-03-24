using System;
using System.Drawing;
using Grasshopper;
using Grasshopper.Kernel;

namespace VoussoirPlugin03.Components.Utils
{
    public class VoussoirPlugin03Info : GH_AssemblyInfo
    {
        public override string Name => "Voussoir";

        //Return a 24x24 pixel bitmap to represent this GHA library.
        public override Bitmap Icon => VoussoirPlugin03.Properties.Resources.Voussoir00;

        //Return a short string describing the purpose of this GHA library.
        public override string Description => "Grasshopper plugin to design stereotomic vaults.";

        public override Guid Id => new Guid("c61bbb89-0359-49eb-a899-bc1f8e69a803");

        //Return a string identifying you or your company.
        public override string AuthorName => "STBIM";

        //Return a string representing your preferred contact details.
        public override string AuthorContact => "";

        //Return a string representing the version.  This returns the same version as the assembly.
        public override string AssemblyVersion => GetType().Assembly.GetName().Version.ToString();
    }

    public class GH_AssemblyPriority
    {
        public static GH_LoadingInstruction PriorityLoad()
        {
            var server = Grasshopper.Instances.ComponentServer;

            server.AddCategoryShortName("Voussoir", "V");
            server.AddCategorySymbolName("Voussoir", 'V');
            server.AddCategoryIcon("Voussoir", VoussoirPlugin03.Properties.Resources.Voussoir00);

            return GH_LoadingInstruction.Proceed;
        }
    }
}