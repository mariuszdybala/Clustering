using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace Cluster
{
    public partial class Form1 : Form
    {
        List<CoordinateBase> Miasta = new List<CoordinateBase>{
            new CoordinateBase(){ Latitude = 50.0467657, Longitude = 20.0048731}, //Kraków
            new CoordinateBase(){ Latitude = 49.9735762, Longitude = 19.8250327}, //Skawina
            new CoordinateBase(){ Latitude = 50.0261438, Longitude = 20.9769061}, //Tarnów
            new CoordinateBase(){ Latitude = 50.4703205, Longitude = 20.7163545}, //Busko Zdrój
            new CoordinateBase(){ Latitude = 50.6395155, Longitude = 20.2893146}, //Jędrzejów
            new CoordinateBase(){ Latitude = 50.2136724, Longitude = 19.0071927}, //Katowice
            new CoordinateBase(){ Latitude = 50.3587496, Longitude = 20.0274582}, //Miechów
            new CoordinateBase(){ Latitude = 50.2857166, Longitude = 19.5617801}}; //Olkusz

        List<CoordinateBase> Clusters;


        public Form1()
        {
            InitializeComponent();
            //GetCenterPoint();
            DrawArea();
        }

        private void DrawArea()
        {
            List<CoordinateBase> UnSortedList = Miasta;
            List<CoordinateBase> SortedList = new List<CoordinateBase>();
            List<CoordinateBase> TempListI = new List<CoordinateBase>();
            List<CoordinateBase> TempListII = new List<CoordinateBase>();
            List<CoordinateBase> TempListIII = new List<CoordinateBase>();
            List<CoordinateBase> TempListIV = new List<CoordinateBase>();

            int UnSortedListCount = UnSortedList.Count;
            var StartPoint = UnSortedList.Min(x => x.Longitude);
            CoordinateBase EndPoint = UnSortedList.First(x => x.Longitude == StartPoint);
            SortedList.Add(UnSortedList.First(x => x.Longitude == StartPoint));
            UnSortedList.Remove(SortedList[0]);

            for (int i = 1; i < UnSortedListCount  ; i++)
            {
                TempListI.Clear();
                TempListII.Clear();
                TempListIII.Clear();
                TempListIV.Clear();
                //szukam punktu najbliższego w I ćw.
                foreach (var miasto in UnSortedList)
                    if (miasto.Latitude > SortedList[i - 1].Latitude && miasto.Longitude > SortedList[i - 1].Longitude)
                        TempListI.Add(miasto);
                    else if (miasto.Latitude < SortedList[i - 1].Latitude && miasto.Longitude > SortedList[i - 1].Longitude)
                        TempListII.Add(miasto);
                    else if (miasto.Latitude < SortedList[i - 1].Latitude && miasto.Longitude < SortedList[i - 1].Longitude)
                        TempListIII.Add(miasto);
                    else if (miasto.Latitude > SortedList[i - 1].Latitude && miasto.Longitude < SortedList[i - 1].Longitude)
                        TempListIV.Add(miasto);


                if (TempListI.Count != 0)
                {
                    StartPoint = TempListI.Min(x => x.Longitude);
                    SortedList.Add(UnSortedList.First(x => x.Longitude == StartPoint));
                    UnSortedList.Remove(SortedList[i]);
                }
                else if (TempListII.Count != 0)
                {
                    StartPoint = TempListII.Max(x => x.Latitude);
                    SortedList.Add(UnSortedList.First(x => x.Latitude == StartPoint));
                    UnSortedList.Remove(SortedList[i]);
                }

                else if (TempListIII.Count != 0)
                {
                    StartPoint = TempListIII.Max(x => x.Longitude);
                    SortedList.Add(UnSortedList.First(x => x.Longitude == StartPoint));
                    UnSortedList.Remove(SortedList[i]);
                }
                else if (TempListIV.Count != 0)
                {
                    StartPoint = TempListIV.Min(x => x.Latitude);
                    SortedList.Add(UnSortedList.First(x => x.Latitude == StartPoint));
                    UnSortedList.Remove(SortedList[i]);
                }
            }
            SortedList.Add(EndPoint);
        }

        private void GetCenterPoint()
        {
            var algo = new AlgorithmClusterBundle();
            List<AlgorithmClusterBundle.LatLong> colorPoints = new List<AlgorithmClusterBundle.LatLong>();
            foreach (var item in Miasta)
                colorPoints.Add(new AlgorithmClusterBundle.LatLong(item.Latitude, item.Longitude, 0,  0));
            List<AlgorithmClusterBundle.Bucket> clusters = algo.RunMarkerClusterer(colorPoints, 50);
            Clusters = new List<CoordinateBase>();
            foreach (var item in clusters)
                Clusters.Add(item.Centroid);

        }
    }
}
