using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Cluster
{
    // By Kunuk Nykjaer
    // Clustering package, Grid clustering, MarkerClusterer and auto increasing K-means clustering
    // should split all class to each file and refactor all global const away
    // quick source code bundle version
    // Version 0.4
    //
    // for visualization this generates a draw.js file which is used by the canvas.html file
    // the canvas.html is distributed elsewhere, e.g. at my blog at http://kunuk.wordpress.com

    public class AlgorithmClusterBundle
    {


        //------------------------------
        // USER CONFIG, CUSTOMIZE BELOW
        public const ClusterType clusterType = ClusterType.MarkerClusterer;

        // use centroid or semi-centroid cluster point placement visualization?
        public const bool DoUpdateAllCentroidsToNearestContainingPoint = false;

        // K-MEANS config
        // heuristic, set tolerance for cluster density, has effect on running time.
        // Set high for many points in dataset, can be lower for fewer points
        //public const double MAX_ERROR = 50;
        public const double MAX_ERROR = 40;

        // MARKERCLUSTERER config
        // box area i.e. cluster size
        //public const int MARKERCLUSTERER_SIZE = 100;
        public static double MARKERCLUSTERER_SIZE = 2;

        // DISTANCECLUSTERER config
        // radius size i.e. cluster size
        //public const int DISTANCECLUSTERER_SIZE = 62;
        private static int DISTANCECLUSTERER_SIZE = 20;

        // USER CONFIG, CUSTOMIZE ABOVE
        //------------------------------



        public enum ClusterType { KMeans, MarkerClusterer, DistanceClusterer } ;

        public static readonly Random Rand = new Random();


        // dictionary lookup key used by grid cluster algo
        public static string GetId(int idx, int idy) //O(1)
        {
            return idx + ";" + idy;
        }


        public List<Bucket> Run(ClusterType clustertype, List<LatLong> points)
        {
            switch (clustertype)
            {
                case ClusterType.KMeans:
                    return new KMeans(points).GetCluster();
                case ClusterType.MarkerClusterer:
                    return new MarkerClusterer(points).GetCluster();
                case ClusterType.DistanceClusterer:
                    return new DistanceClusterer(points).GetCluster();
            }
            return null;
        }

        public List<Bucket> RunMarkerClusterer(List<LatLong> points, int bufferKm)
        {
            MARKERCLUSTERER_SIZE = (double)bufferKm / 50;
            return new MarkerClusterer(points).GetCluster();
        }

        public List<Bucket> RunDistanceClusterer(List<LatLong> points, int bufferKm)
        {
            DISTANCECLUSTERER_SIZE = bufferKm;
            return new DistanceClusterer(points).GetCluster();
        }




        [Serializable()]
        public class LatLong : CoordinateBase, IComparable
        {

            public int Size { get; set; }
            public double Numerator { get; set; }
            public double Denominator { get; set; }

            public LatLong() { }
            public LatLong(double latitude, double longitude)
            {
                Latitude = latitude;
                Longitude = longitude;
            }

            public LatLong(double latitude, double longitude, double numerator, double denominator)
            {
                Latitude = latitude;
                Longitude = longitude;
                Numerator = numerator;
                Denominator = denominator;
            }

            public LatLong(LatLong p) //clone
            {
                this.Latitude = p.Latitude;
                this.Longitude = p.Longitude;
                this.Size = p.Size;
                this.Numerator = p.Numerator;
                this.Denominator = p.Denominator;
            }

            public int CompareTo(object o) // if used in sorted list
            {
                if (this.Equals(o))
                    return 0;

                var other = (LatLong)o;
                if (this.Latitude > other.Latitude)
                    return -1;
                if (this.Latitude < other.Latitude)
                    return 1;

                return 0;
            }

            // used by k-means random distinct selection of cluster point
            public override int GetHashCode()
            {
                var x = Latitude * 10000; //make the decimals be important
                var y = Longitude * 10000;
                var r = x * 17 + y * 37;
                return (int)r;
            }
            private const int ROUND = 6;
            public override bool Equals(Object o)
            {
                if (o == null)
                    return false;
                var other = o as LatLong;
                if (other == null)
                    return false;

                // rounding could be skipped
                // depends on granularity of wanted decimal precision
                // note, 2 points with same x,y is regarded as being equal
                var x = Math.Round(this.Latitude, ROUND) == Math.Round(other.Latitude, ROUND);
                var y = Math.Round(this.Longitude, ROUND) == Math.Round(other.Longitude, ROUND);
                return x && y;
            }

        }
        public class Bucket
        {
            public string Id { get; private set; }
            public List<LatLong> Points { get; private set; }
            public LatLong Centroid { get; set; }
            public int Idx { get; private set; }
            public int Idy { get; private set; }
            public double ErrorLevel { get; set; } // clusterpoint and points avg dist
            private bool _IsUsed;
            public bool IsUsed
            {
                get { return _IsUsed && Centroid != null; }
                set { _IsUsed = value; }
            }
            public Bucket(string id)
            {
                IsUsed = true;
                Centroid = null;
                Points = new List<LatLong>();
                Id = id;
            }
            public Bucket(int idx, int idy)
            {
                IsUsed = true;
                Centroid = null;
                Points = new List<LatLong>();
                Idx = idx;
                Idy = idy;
                Id = GetId(idx, idy);
            }
        }


        public abstract class BaseClusterAlgorithm
        {
            public List<LatLong> BaseDataset; // all points
            //id, bucket
            public readonly Dictionary<string, Bucket> BaseBucketsLookup =
                new Dictionary<string, Bucket>();

            public BaseClusterAlgorithm() { }
            public BaseClusterAlgorithm(List<LatLong> dataset)
            {
                if (dataset == null || dataset.Count == 0)
                    throw new ApplicationException(
                        string.Format("dataset is null or empty"));

                BaseDataset = dataset;
            }

            public abstract List<Bucket> GetCluster();
            //O(k?? random fn can be slow, but is not slow because atm the k is always 1)
            public static LatLong[] BaseGetRandomCentroids(List<LatLong> list, int k)
            {
                var set = new HashSet<LatLong>();
                int i = 0;
                var kcentroids = new LatLong[k];

                int MAX = list.Count;
                while (MAX >= k)
                {
                    int index = Rand.Next(0, MAX - 1);
                    var xy = list[index];
                    if (set.Contains(xy))
                        continue;

                    set.Add(xy);
                    kcentroids[i++] = new LatLong(xy.Latitude, xy.Longitude);

                    if (i >= k)
                        break;
                }
                return kcentroids;
            }

            public List<Bucket> BaseGetClusterResult()
            {

                // collect used buckets and return the result
                var clusterPoints = new List<Bucket>();
                foreach (var item in BaseBucketsLookup)
                {
                    var bucket = item.Value;
                    if (bucket.IsUsed)
                    {
                        bucket.Centroid.Size = bucket.Points.Count;
                        clusterPoints.Add(bucket);
                    }
                }

                return clusterPoints;
            }
            public static LatLong BaseGetCentroidFromCluster(List<LatLong> list) //O(n)
            {
                int count = list.Count;
                if (list == null || count == 0)
                    return null;

                // color is set for the points and the cluster point here
                LatLong centroid = new LatLong(0, 0) { Size = list.Count };//O(1)
                foreach (LatLong p in list)
                {
                    centroid.Latitude += p.Latitude;
                    centroid.Longitude += p.Longitude;
                }
                centroid.Latitude /= count;
                centroid.Longitude /= count;
                var cp = new LatLong(centroid.Latitude, centroid.Longitude) { Size = count };

                return cp;
            }
            //O(k*n)
            public static void BaseSetCentroidForAllBuckets(IEnumerable<Bucket> buckets)
            {
                foreach (var item in buckets)
                {
                    var bucketPoints = item.Points;
                    var cp = BaseGetCentroidFromCluster(bucketPoints);
                    item.Centroid = cp;
                }
            }
            public double BaseGetTotalError()//O(k)
            {
                int centroidsUsed = BaseBucketsLookup.Values.Count(b => b.IsUsed);
                double sum = BaseBucketsLookup.Values.
                    Where(b => b.IsUsed).Sum(b => b.ErrorLevel);
                return sum / centroidsUsed;
            }
            public string BaseGetMaxError() //O(k)
            {
                double maxError = -double.MaxValue;
                string id = string.Empty;
                foreach (var b in BaseBucketsLookup.Values)
                {
                    if (!b.IsUsed || b.ErrorLevel <= maxError)
                        continue;

                    maxError = b.ErrorLevel;
                    id = b.Id;
                }
                return id;
            }
            public LatLong BaseGetClosestPoint(LatLong from, List<LatLong> list) //O(n)
            {
                double min = double.MaxValue;
                LatLong closests = null;
                foreach (var p in list)
                {
                    var d = MathTool.Distance(from, p);
                    if (d >= min)
                        continue;

                    // update
                    min = d;
                    closests = p;
                }
                return closests;
            }
            public LatLong BaseGetLongestPoint(LatLong from, List<LatLong> list) //O(n)
            {
                double max = -double.MaxValue;
                LatLong longest = null;
                foreach (var p in list)
                {
                    var d = MathTool.Distance(from, p);
                    if (d <= max)
                        continue;

                    // update
                    max = d;
                    longest = p;
                }
                return longest;
            }
            // assign all points to nearest cluster
            public void BaseUpdatePointsByCentroid()//O(n*k)
            {
                int count = BaseBucketsLookup.Count();

                // clear points in the buckets, they will be re-inserted
                foreach (var bucket in BaseBucketsLookup.Values)
                    bucket.Points.Clear();

                foreach (LatLong p in BaseDataset)
                {
                    double minDist = Double.MaxValue;
                    string index = string.Empty;
                    //for (int i = 0; i < count; i++)
                    foreach (var i in BaseBucketsLookup.Keys)
                    {
                        var bucket = BaseBucketsLookup[i];
                        if (bucket.IsUsed == false)
                            continue;

                        var centroid = bucket.Centroid;
                        var dist = MathTool.Distance(p, centroid);
                        if (dist < minDist)
                        {
                            // update
                            minDist = dist;
                            index = i;
                        }
                    }
                    //update color for point to match centroid and re-insert
                    var closestBucket = BaseBucketsLookup[index];
                    closestBucket.Points.Add(p);
                }
            }

            // update centroid location to nearest point, 
            // e.g. if you want to show cluster point on a real existing point area
            //O(n)
            public void BaseUpdateCentroidToNearestContainingPoint(Bucket bucket)
            {
                if (bucket == null || bucket.Centroid == null ||
                    bucket.Points == null || bucket.Points.Count == 0)
                    return;

                var closest = BaseGetClosestPoint(bucket.Centroid, bucket.Points);
                bucket.Centroid.Latitude = closest.Latitude;
                bucket.Centroid.Longitude = closest.Longitude;
            }
            //O(k*n)
            public void BaseUpdateAllCentroidsToNearestContainingPoint()
            {
                foreach (var bucket in BaseBucketsLookup.Values)
                    BaseUpdateCentroidToNearestContainingPoint(bucket);
            }
        }

        public class DistanceClusterer : BaseClusterAlgorithm
        {
            public DistanceClusterer(List<LatLong> dataset)
                : base(dataset)
            {
            }

            public override List<Bucket> GetCluster()
            {
                var cluster = RunClusterAlgo();
                return cluster;
            }

            // O(k*n)
            List<Bucket> RunClusterAlgo()
            {
                // put points in buckets     
                int allPointsCount = BaseDataset.Count;
                var firstPoint = BaseDataset[0];
                var firstId = 0.ToString();
                var firstBucket = new Bucket(firstId) { Centroid = firstPoint };
                BaseBucketsLookup.Add(firstId, firstBucket);

                for (int i = 1; i < allPointsCount; i++)
                {
                    var set = new HashSet<string>(); //cluster candidate list
                    var p = BaseDataset[i];
                    // iterate clusters and collect candidates
                    foreach (var bucket in BaseBucketsLookup.Values)
                    {
                        var isInCluster = MathTool.DistWithin(p, bucket.Centroid, DISTANCECLUSTERER_SIZE);
                        if (!isInCluster)
                            continue;

                        set.Add(bucket.Id);
                        //use first, short dist will be calc at last step before returning data
                        break;
                    }

                    // if not within box area, then make new cluster   
                    if (set.Count == 0)
                    {
                        var pid = i.ToString();
                        var newbucket = new Bucket(pid) { Centroid = p };
                        BaseBucketsLookup.Add(pid, newbucket);
                    }
                }

                //important, align all points to closest cluster point
                BaseUpdatePointsByCentroid();

                return BaseGetClusterResult();
            }

        }

        public class MarkerClusterer : BaseClusterAlgorithm
        {
            public MarkerClusterer(List<LatLong> dataset)
                : base(dataset)
            {
            }

            public override List<Bucket> GetCluster()
            {
                var cluster = RunClusterAlgo();
                return cluster;
            }

            // O(k*n)
            List<Bucket> RunClusterAlgo()
            {
                // put points in buckets     
                int allPointsCount = BaseDataset.Count;
                var firstPoint = BaseDataset[0];
                var firstId = 0.ToString();
                var firstBucket = new Bucket(firstId) { Centroid = firstPoint };
                BaseBucketsLookup.Add(firstId, firstBucket);

                for (int i = 1; i < allPointsCount; i++)
                {
                    var set = new HashSet<string>(); //cluster candidate list
                    var p = BaseDataset[i];
                    // iterate clusters and collect candidates
                    foreach (var bucket in BaseBucketsLookup.Values)
                    {
                        var isInCluster = MathTool.BoxWithin(p, bucket.Centroid, MARKERCLUSTERER_SIZE);
                        if (!isInCluster)
                            continue;

                        set.Add(bucket.Id);
                        //use first, short dist will be calc at last step before returning data
                        break;
                    }

                    // if not within box area, then make new cluster   
                    if (set.Count == 0)
                    {
                        var pid = i.ToString();
                        var newbucket = new Bucket(pid) { Centroid = p };
                        BaseBucketsLookup.Add(pid, newbucket);
                    }
                }

                //important, align all points to closest cluster point
                BaseUpdatePointsByCentroid();

                return BaseGetClusterResult();
            }
        }

        // O(exponential) ~ can be slow when n or k is big
        public class KMeans : BaseClusterAlgorithm
        {
            private readonly int _InitClusterSize; // start from this cluster points
            // Rule of thumb k = sqrt(n/2)        

            // cluster point optimization iterations
            private const int _MaxIterations = 100;
            private const int _MaxClusters = 100;

            public KMeans(List<LatLong> dataset)
                : base(dataset)
            {
                _InitClusterSize = 1;
            }

            public override List<Bucket> GetCluster()
            {
                var cluster = RunClusterAlgo();
                return cluster;
            }

            List<Bucket> RunClusterAlgo()
            {
                /*
                ITERATIVE LINEAR ADDING CLUSTER UNTIL 
                REPEATE iteration of clusters convergence 
                   until the max error is small enough
                if not, insert a new cluster at worst place 
                (farthest in region of worst cluster area) and run again 
                keeping current cluster points              
             
                 // one iteration of clusters convergence is defined as ..
              1) Random k centroids
              2) Cluster data by euclidean distance to centroids
              3) Update centroids by clustered data,     
              4) Update cluster
              5) Continue last two steps until error is small, error is sum of diff
                  between current and updated centroid                          
               */

                RunAlgo();

                if (DoUpdateAllCentroidsToNearestContainingPoint)
                    BaseUpdateAllCentroidsToNearestContainingPoint();
                return BaseGetClusterResult();
            }

            void RunAlgo()
            {
                // Init clusters
                var centroids = BaseGetRandomCentroids(BaseDataset, _InitClusterSize);
                for (int i = 0; i < centroids.Length; i++)
                {
                    var pid = i.ToString();
                    var newbucket = new Bucket(pid) { Centroid = centroids[i] };
                    BaseBucketsLookup.Add(pid, newbucket);
                }

                //
                double currentMaxError = double.MaxValue;
                while (currentMaxError > MAX_ERROR && BaseBucketsLookup.Count < _MaxClusters)
                {
                    RunIterationsUntilKClusterPlacementAreDone();

                    var id = BaseGetMaxError();
                    var bucket = BaseBucketsLookup[id];
                    currentMaxError = bucket.ErrorLevel; //update
                    if (currentMaxError > MAX_ERROR)
                    {
                        // Here it is linear speed when putting one new centroid at a time
                        // should be semi-fast because the new point is inserted at best area 
                        //from current centroids view.
                        // possible improvement exists by log2 search by inserting multiple centroids and 
                        // reducing centroid again if needed

                        // put new centroid in area where maxError but farthest away from current centroid in area
                        var longest = BaseGetLongestPoint(bucket.Centroid, bucket.Points);
                        var newcentroid = new LatLong(longest);
                        var newid = BaseBucketsLookup.Count.ToString();
                        var newbucket = new Bucket(newid) { Centroid = newcentroid };
                        BaseBucketsLookup.Add(newid, newbucket);
                    }
                }
            }

            void RunIterationsUntilKClusterPlacementAreDone()
            {
                double prevError = Double.MaxValue;
                double currError = Double.MaxValue;

                for (int i = 0; i < _MaxIterations; i++)
                {
                    prevError = currError;
                    currError = RunOneIteration();
                    if (currError >= prevError) // no improvement
                        break;
                }
            }

            double RunOneIteration() //O(k*n)
            {

                // update points, assign points to cluster
                BaseUpdatePointsByCentroid();
                // update centroid pos by its points
                BaseSetCentroidForAllBuckets(BaseBucketsLookup.Values);//O(k*n)
                var clustersCount = BaseBucketsLookup.Count;
                for (int i = 0; i < clustersCount; i++)
                {
                    var currentBucket = BaseBucketsLookup[i.ToString()];
                    if (currentBucket.IsUsed == false)
                        continue;

                    //update centroid                
                    var newcontroid = BaseGetCentroidFromCluster(currentBucket.Points);
                    //no need to update color, autoset
                    currentBucket.Centroid = newcontroid;
                    currentBucket.ErrorLevel = 0;
                    //update error                
                    foreach (var p in currentBucket.Points)
                    {
                        var dist = MathTool.Distance(newcontroid, p);
                        currentBucket.ErrorLevel += dist;
                    }
                    var val = currentBucket.ErrorLevel / currentBucket.Points.Count;
                    currentBucket.ErrorLevel = val; //Math.Sqrt(val);
                }

                return BaseGetTotalError();
            }

        }

        public class MathTool
        {

            private const double Exp = 2; // 2=euclid, 1=manhatten

            public static double Distance(LatLong a, LatLong b)
            {
                return StraightLineDistance(a.Latitude, a.Longitude, b.Latitude, b.Longitude);
                //return Math.Pow(Math.Pow(Math.Abs(a.X - b.X), Exp) +
                //    Math.Pow(Math.Abs(a.Y - b.Y), Exp), 1.0 / Exp);
            }


            public static double Min(double a, double b)
            {
                return a <= b ? a : b;
            }
            public static double Max(double a, double b)
            {
                return a >= b ? a : b;
            }

            public static bool DistWithin(LatLong a, LatLong b, double d)
            {
                // var dist = Distance(a, b);
                var dist = StraightLineDistance(a.Latitude, a.Longitude, b.Latitude, b.Longitude);
                return dist < d;
            }

            /// <summary>
            /// Calculates the distance between two coordinates
            /// </summary>
            /// <param name="Lat1">Latitude of first coordinate</param>
            /// <param name="Long1">Longitude of first coordinate</param>
            /// <param name="Lat2">Latitude of second coordinate</param>
            /// <param name="Long2">Longitude of second coordinate</param>
            /// <returns>Number of km</returns>
            public static double StraightLineDistance(double Lat1, double Long1, double Lat2, double Long2)
            {
                double distance = 0;
                double x = 0;
                double y = 0;

                x = 69.1 * (Lat1 - Lat2);
                y = 69.1 * (Long1 - Long2) * System.Math.Cos(Lat2 / 57.3);

                //calculation base : Miles
                distance = System.Math.Sqrt(x * x + y * y);

                //Distance calculated in Kilometres
                return distance * 1.609;

                // http://www.codeproject.com/KB/cs/distancebetweenlocations.aspx
                /*
                    The Haversine formula according to Dr. Math.
                    http://mathforum.org/library/drmath/view/51879.html
                
                    dlon = lon2 - lon1
                    dlat = lat2 - lat1
                    a = (sin(dlat/2))^2 + cos(lat1) * cos(lat2) * (sin(dlon/2))^2
                    c = 2 * atan2(sqrt(a), sqrt(1-a)) 
                    d = R * c
                
                    Where
                        * dlon is the change in longitude
                        * dlat is the change in latitude
                        * c is the great circle distance in Radians.
                        * R is the radius of a spherical Earth.
                        * The locations of the two points in 
                            spherical coordinates (longitude and 
                            latitude) are lon1,lat1 and lon2, lat2.
                
                double dDistance = Double.MinValue;
                double dLat1InRad = Lat1 * (Math.PI / 180.0);
                double dLong1InRad = Long1 * (Math.PI / 180.0);
                double dLat2InRad = Lat2 * (Math.PI / 180.0);
                double dLong2InRad = Long2 * (Math.PI / 180.0);

                double dLongitude = dLong2InRad - dLong1InRad;
                double dLatitude = dLat2InRad - dLat1InRad;

                // Intermediate result a.

                double a = Math.Pow(Math.Sin(dLatitude / 2.0), 2.0) +
                           Math.Cos(dLat1InRad) * Math.Cos(dLat2InRad) *
                           Math.Pow(Math.Sin(dLongitude / 2.0), 2.0);

                // Intermediate result c (great circle distance in Radians).

                double c = 2.0 * Math.Atan2(Math.Sqrt(a), Math.Sqrt(1.0 - a));


                // const Double kEarthRadiusMiles = 3956.0;

                const Double kEarthRadiusKms = 6376.5;
                dDistance = kEarthRadiusKms * c;
                //return (int)Math.Round(dDistance, 0);
                return dDistance;
                 */
            }

            public static bool BoxWithin(LatLong a, LatLong b, double boxsize)
            {
                var d = boxsize / 2;
                var withinX = a.Latitude - d <= b.Latitude && a.Latitude + d >= b.Latitude;
                var withinY = a.Longitude - d <= b.Longitude && a.Longitude + d >= b.Longitude;
                return withinX && withinY;
            }
        }

    }
}
