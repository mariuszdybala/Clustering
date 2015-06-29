using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Cluster
{
  public  class CoordinateBase
    {
        public double Latitude { get; set; }
        public double Longitude { get; set; }
        public int Potencial { get; set; }
        public int Size { get; set; }
        private static readonly Random getrandom = new Random();
        private static readonly object syncLock = new object();

      public CoordinateBase()
        {
            Potencial = GetRandomNumber(1, 100);
        }
        //public void SetSizeBubble()
        //{
        //    Size = MyFunc();
        //}

      public static int GetRandomNumber(int min, int max)
      {
          lock (syncLock)
          { // synchronize
              return getrandom.Next(min, max);
          }
      }
    }
}
