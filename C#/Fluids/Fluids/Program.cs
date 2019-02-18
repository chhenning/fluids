using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Fluids
{
    class Program
    {
        static void Main(string[] args)
        {
            Stopwatch stopWatch = new Stopwatch();
            stopWatch.Start();


            Cube c = new Cube(1024, 0, 0, 0.1f);

            c.AddDensity(100, 100, 100);
            c.Step();

            stopWatch.Stop();
            TimeSpan ts = stopWatch.Elapsed;

            Console.WriteLine(string.Format("{0}", stopWatch.ElapsedMilliseconds));
        }
    }
}
