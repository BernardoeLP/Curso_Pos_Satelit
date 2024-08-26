using System.Globalization;
using Numpy;
//using NumpyDotNet;

string nl = Environment.NewLine;
double c = 299792458;                           // m/s  de ITRF
double μ = 3.986005E14;                         // m3/s2  Earth gravitational constant
double ωe = 7.2921151467E-5;                    // radians/s Angular Velocity of the Earth
double[] ωE = new double[3] { 0.0, 0.0, ωe };   // Velocity of the Earth (vector)

DateTime tGPS0 = new(1980, 1, 6, 0, 0, 0);
// el momento en que se realiza la observación
DateTime tobs = new(2001, 3, 19, 0, 15, 0);
Int32 tobs_seg = Convert.ToInt32((tobs - tGPS0).TotalSeconds);
/*            
List<double> L = new List<double>();
List<double[]> A = new List<double[]>();

var L = np.empty((0, 1));
var A = np.empty((0, 2));
*/

NDarray L = np.empty((0, 1));
NDarray A = np.empty((0, 2));

//static double[] X1 = new double[4] { 0.0, 0.0, 0.0, 0.0 };
NDarray X1 = np.array(new double[4]);
List<double> C = new List<double>();

List<string> satord = new List<string>();  // es un 'parche' sólo para mostrar a que satélite cooresponde cada correccion

Dictionary<string, double[]> P2D           // Observables incluyendo código P
        = new Dictionary<string, double[]>()
{                       //             C/A          Fase L1           Fase L2           P1           P2
            { "28", new double[]{ 23334619.807 ,-10956241.60549 , -8525575.42946 ,23334619.277 ,23334623.308} },
            { "13", new double[]{ 22586189.129 ,-12029006.00949 , -9358589.69746 ,22586189.572 ,22586192.629} },
            { "01", new double[]{ 25167667.280 ,   -19838.71849 ,   275138.51344 ,25167667.160 ,25167670.651} },
            { "27", new double[]{ 20873169.266 ,-23420787.65349 ,-18223955.59247 ,20873168.733 ,20873173.057} },
            { "24", new double[]{ 23371141.291 , -8542099.54349 , -6138557.63846 ,23371140.735 ,23371147.385} },
            { "10", new double[]{ 21505922.486 ,-16684416.38149 ,-12470663.29047 ,21505921.442 ,21505925.535} },
            { "08", new double[]{ 20958408.428 ,-22772002.73249 ,-17730879.20747 ,20958407.438 ,20958412.414} }
};


double f1 = 1575.42;
double f12 = f1 * f1;
double f2 = 1227.60;
double f22 = f2 * f2;
double fdif = f12 - f22;
Dictionary<string, double> PD = new Dictionary<string, double>();
Dictionary<string, double> TT = new Dictionary<string, double>();

foreach (var key in P2D.Keys)
{
    string s = key.ToString();
    double PseudoD = (f12 * P2D[s][3] - f22 * P2D[s][4]) / fdif;
    PD[s] = PseudoD;
    // Obtengo tiempo de tránsito [s]
    TT[s] = PseudoD / c;
    //Console.WriteLine("Sat {0}   PD: {1}    TT: {2}", s, PD[s], TT[s]);
}

int cant_sat = PD.Count();

Dictionary<string, double[]> Precisas
    = new Dictionary<string, double[]>()        // Efemérides Precisas [Km] , [us]
    {
                                //        X             Y              Z           clk
                { "01" , new double[]{  581.886423, 25616.666528 ,  7088.545471, 169.092800} },
                { "08" , new double[]{22018.953984,  2878.718252 , 14451.124018,   9.709180} },
                { "10" , new double[]{10103.948910, -10925.429662, 22009.912003,   1.148951} },
                { "13" , new double[]{ 7525.432597,  20488.591201, 15216.097471,  -0.655216} },
                { "24" , new double[]{22368.646126, -12657.086060,  6934.928617,  36.698468} },
                { "27" , new double[]{15057.427636,   9402.947329, 20171.667340,  14.763242} },
                { "28" , new double[]{-5895.039751,  14576.928529, 21538.074040,  14.267922} }
    };

Dictionary<string, double[]> Calculadas = new();

Dictionary<string, List<Efemerides>> mensajes = new Dictionary<string, List<Efemerides>>();
Dictionary<string, DateTime> tefem = new Dictionary<string, DateTime>();


List<string> satlist = new List<string>();

Dictionary<string, double> Delta_ant
    = new Dictionary<string, double>()         // Corrección de centro de fase de antena
    {
                            { "28" ,  1.04280},                     //     para cada sat [m]
                            { "13" ,  1.38950},
                            { "01" ,  2.38080},
                            { "27" ,  2.63340},
                            { "24" ,  2.60380},
                            { "10" ,  2.54650},
                            { "08" ,  2.57810}
    };

double[] Estacion = new double[4] {            // Coord. precisas de la Estación [m]
                            3370658.6942,     // X
                            711877.0150,     // Y
                            5349786.8637,     // Z
                                    0.0         // c * Delta t
                            };
double[] Coord = new  double[4];               // Coord. calculadas de la Estación [m]
double[] rr = new double[3];

static double newton(Func<double, double> fn, Func<double, double> Df, double x0, double epsilon, int max_iter)
{
    /*Approximate solution of f(x)=0 by Newton's method.

    Parameters
    ----------
    fn: function
        Function for which we are searching for a solution f(x) = 0.
    Df : function
        Derivative of f(x).
    x0 : number
        Initial guess for a solution f(x) = 0.
    epsilon : number
        Stopping criteria is abs(f(x)) < epsilon.
    max_iter : integer
        Maximum number of iterations of Newton's method.

    Returns
    ------ -
    xn : number
        Implement Newton's method: compute the linear approximation
        of f(x) at xn and find x intercept by the formula
            x = xn - f(xn) / Df(xn)
        Continue until abs(f(xn)) < epsilon and return xn.
        If Df(xn) == 0, return None.If the number of iterations
        exceeds max_iter, then return None.
            Examples
    --------
    >>> f = lambda x: x * *2 - x - 1
    >>> Df = lambda x: 2 * x - 1
    >>> newton(f, Df, 1, 1e-8, 10)
    Found solution after 5 iterations.
    1.618033988749989
    */
    double xn = x0;
    for (int j = 0; j < max_iter; j++)
    {
        double fxn = fn(xn);
        if (Math.Abs(fxn) < epsilon)
            //Console.WriteLine('Found solution after',n,'iterations.')
            return xn;
        double Dfxn = Df(xn);
        if (Dfxn == 0)
        {
            Console.WriteLine("Zero derivative. No solution found.");
            return Double.NaN;
        };
        xn = xn - fxn / Dfxn;
    };
    Console.WriteLine("Exceeded maximum iterations. No solution found.");
    return Double.NaN;
}

Dictionary<string, double[]> coord_obs()
{
    // Arma la matriz con coordenadas de cada satélite
    //    en el tiempo de observación

    Dictionary<string, double[]> Calc = new();
    foreach (var sati in TT.Keys)
    {
        List <Efemerides> s = mensajes[sati];
        foreach (Efemerides h in s)
        {
            if (h.FechaHora == tefem[sati])
            {
                double tefem_seg = ((tefem[sati] - tGPS0).TotalSeconds);
                double a = h.sqrtA * h.sqrtA;
                double a_cubo = a * a * a;

                // tiempo entre la efemeride transmitida y
                // el momento en que quiero saber las coordenadas del satélite
                double Delta_t = (tobs_seg - TT[sati]) - tefem_seg;


                double M = h.M0 + (Math.Sqrt(μ / a_cubo) + h.Delta_n) * Delta_t;

                Func<double, double> fE = x => M + h.e * Math.Sin(x) - x;
                Func<double, double> dfE = x => h.e * Math.Cos(x) - 1;
                double E = newton(fE, dfE, 0, 1E-6, 4);
                Console.WriteLine("Satelite: "+sati+string.Format("  E: {0,15:F10}",E));

                double r0 = a * (1 - h.e * Math.Cos(E));
                double f = 2 * Math.Atan(Math.Sqrt(1 + h.e) / Math.Sqrt(1 - h.e) * Math.Tan(E / 2));
                double u0 = h.omega + f;


                double Ω = h.OMEGA + h.OMEGA_DOT * Delta_t;
                double ω = h.omega + h.Cuc * Math.Cos(2 * u0) + h.Cus * Math.Sin(2 * u0);
                double r = r0 + h.Crc * Math.Cos(2 * u0) + h.Crs * Math.Sin(2 * u0);
                double i = h.i0 + h.Cic * Math.Cos(2 * u0) + h.Cis * Math.Sin(2 * u0) + h.idot * Delta_t;
                double ϴ = ωe * (tobs_seg - TT[sati]);
                double u = ω + f;
                r += Delta_ant[sati];

                double[,] xyz = new double[3, 3]
                    {   {r * Math.Cos(u), 0, 0},
                        {r * Math.Sin(u), 0, 0},
                        {            0  , 0, 0}};

                double[,] R31 = new double[3, 3]
                    {   {  Math.Cos(-Ω), Math.Sin(-Ω), 0},
                        { -Math.Sin(-Ω), Math.Cos(-Ω), 0},
                        {          0   ,         0   , 1}};

                double[,] R32 = new double[3, 3]
                    {   {  Math.Cos(ϴ), Math.Sin(ϴ),  0},
                        { -Math.Sin(ϴ), Math.Cos(ϴ),  0},
                        {          0  ,         0  ,  1}};

                double[,] R1 = new double[3, 3]
                    {   { 1,         0   ,        0    },
                        { 0, Math.Cos(-i), Math.Sin(-i)},
                        { 0,-Math.Sin(-i), Math.Cos(-i)}};

                NDarray R3R3 =  np.matmul( np.array(R32),  np.array(R31));
                NDarray RR   =  np.matmul(R3R3,  np.array(R1));


                NDarray xyz_prima = np.matmul(RR,  np.array(xyz));
                double x = (double)xyz_prima[0][0];
                double y = (double)xyz_prima[1][0];
                double z = (double)xyz_prima[2][0];
                Calc[sati] = new double[4]{ x, y, z, Precisas[sati][3]};
            }
        }
    }
    return Calc;
}

void arma_matriz(ref NDarray a, ref NDarray l)
{
    /*
    """     al armar la matriz de diseño se puede poner c *Delta_t en la
            4ta.columna en vez de C, así los resultados
            dan en Distancia, en lugar de Delta_t.
            Para eso, pongo todos unos en vez de c,
            así la incógnita incluye a c y quedan números más manejables
    """
    */
    C.Clear();
    List<double> ll = new List<double>();
    satord.Clear();
    //int j = 0;
    bool primera = true;
    foreach (var st in Calculadas.Keys)
    {
        double[] s = Calculadas[st];
        double[] rs = s[0..3];       // Coordenadas (x, Y, Z) calculadas del satélite "st"
        double dX = Coord[0] - s[0];
        double dY = Coord[1] - s[1];
        double dZ = Coord[2] - s[2];
        double ρ = Math.Sqrt(dX * dX + dY * dY + dZ * dZ);
        double ρSagnac = (double) np.dot(np.subtract(rr, rs), np.cross(ωE, rr)) / c;

        Console.WriteLine(string.Format("Sat: {0:s}  Coor. x Sagnac: {1,10:F5}", st, ρSagnac));

        double[] fila = new double[4]{ (dX / ρ), (dY / ρ), (dZ / ρ), 1.0 };  // opcion con incognita c * Delta_t

        NDarray d = np.array(fila);
        if (primera)
        {
            a = np.append(a, d);
            primera = false;
        }
        else
            a = np.vstack( a, np.array(fila));

        // diferencia Observado - Calculado
        double err = PD[st] - ρ - Coord[3] + c * s[3] / 1E6 - ρSagnac;
        // Si s[3] > 0 el satélite atrasa con respecto a GPS time, entonces "sumo" error en distancia

        //l = np.append(l,np.array(err));
        ll.Add(err);
        satord.Add(st);
        /*
        linea_C =[0 for i in range(cant_sat)];
        linea_C[j] = err * err;
        j += 1;
        C.append(linea_C);
        */
    }
    l = np.array(ll.ToArray());
}

void imprime_dif_sat()
{
    foreach (var s in Calculadas.Keys)
    {
        Console.WriteLine(nl + "Sat: " + s);
        string linea1 = " Calculadas : ";
        string linea2 = " Precisas   : ";
        string linea3 = " Diferencia : ";
        double d;
        double mod = 0.0;
        for (int j=0; j< 3; j++)
        {
            linea1 += string.Format("  {0,15:F5} m",Calculadas[s][j]);
            linea2 += string.Format("  {0,15:F5} m",Precisas[s][j] * 1000);
            d = Calculadas[s][j] - Precisas[s][j] * 1000;
            linea3 += string.Format("  {0,15:F5} m",d);
            mod += d * d;
        }
        Console.WriteLine(linea1);
        Console.WriteLine(linea2);
        Console.WriteLine(linea3);
        Console.WriteLine(string.Format("  Módulo de la diferencia: {0,15:F5} m",Math.Sqrt(mod)));
    }
}

void imprime_resu()
{
    //  Imprime resultados parciales

    Console.WriteLine(nl + "Observado-calculado");
    Console.WriteLine("SAT         Dif" + nl);
    for (int j = 0; j < L.size; j++)
        Console.WriteLine(satord[j] + string.Format("   {0,15:F6}  ", L[j]));
    Console.WriteLine();
    Console.WriteLine(string.Format("Dispersión:  {0,8:F6}",np.std(L))+ nl);

    /*
    for dif in L:
        Console.WriteLine("{:15.6f}  ".format(dif))
    Console.WriteLine()

    for linea in P:
        Console.WriteLine(linea)
    Console.WriteLine()
    */
    
    Console.WriteLine(nl+"Matriz de diseño"+nl);
    for (int i = 0; i < cant_sat; i++)
    {
        string reng = "";
        for (int j=0;j<A[i].size;j++)
                reng += string.Format("{0,20:F16}  ",A[i][j]);
        Console.WriteLine(reng);
    }
    
    //Console.WriteLine("\n")
    /*
    Console.WriteLine()
    Console.WriteLine("Coord Calculada")
    linea = ""
    for i in Coord:
        linea += "{:20.16f}  ".format(i)
    Console.WriteLine(linea)
    Console.WriteLine()
    */
    Console.WriteLine(" - - - - - - - - - -");
    Console.WriteLine("Diferencias Calculadas para corregir las coordenadas");
    Console.WriteLine("           X                   Y                   Z                   t");
    string linea = "";
    for (int i = 0; i < X1.size;i++)
        linea += string.Format("{0,20:F16}  ", X1[i]);
    Console.WriteLine(linea);
    NDarray dift = np.empty((3,1));

    // Copying the contents of 'x' to array 'y'
    for (int i=0;i<3;i++)
        dift[i] = X1[i];
    //NDarray dift =  X1[0..3];
    Console.WriteLine(nl + string.Format("Módulo de la corrección calculada: {0:20.16f} m", np.linalg.norm(dift)));
    NDarray diftr = X1;
    diftr[3] = X1[3] * c;
    Console.WriteLine("Módulo de la corrección calculada" + nl + string.Format("              incluyendo el reloj: {0,10:F6} m", np.linalg.norm(diftr)) + nl);
    Console.WriteLine(" - - - - - - - - - -" + nl);
}

void imprime_Correg()
{
    // Imprime una vez recalculado
    Console.WriteLine("Coord Corregida,    Precisa,                Diferencia (Calculada - Precisa)");
    string linea = "";
    int j = 0;
    double acu = 0;
    foreach(double i in Coord)
        {
        if (j < 3)
        {
            linea = string.Format("{0,15:F5}  {1,15:F5}  -->{2,15:F5} m",i, Estacion[j], i - Estacion[j]);
            acu += (i - Estacion[j]) * (i - Estacion[j]);
        }
        else
        {
            linea = nl + string.Format(" Dif. entre sitios: {0,15:F5} m", Math.Sqrt(acu)) + nl;
            linea += string.Format(" Delta_t:           {0,15:F5} useg",i / c * 1E6);
            linea += nl + string.Format(" c * Delta_t:       {0,15:F5} m",i) + nl;
        }
        j++ ;
        Console.WriteLine(linea);
        }
}

double calcula_angulo_dif()
{
    double[] vect_dif = new double[3];
    double[] Coord_est = Estacion[0..3];
    for (int j = 0; j < 3; j++)
        vect_dif[j] = Coord[j] - Estacion[j];
    double modprod = (double) (np.linalg.norm(Coord_est) * np.linalg.norm(vect_dif));
    return Math.Acos((double) np.dot(Coord_est, vect_dif) / modprod) * 180 / Math.PI ;
}
    
//  / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / -
// Para tomar como coordenadas a-priori de la estacion.
//    se le agrega una diferencia a las precisas
Random random = new Random();
for (int i = 0; i < 3; i++) Coord[i] = Estacion[i] + (random.NextDouble() - 0.5) * 5000;
Coord[3] =0.0;
// ------------------------------------------------------------------------------------
// # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
//  / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / -
// Levanto las efemérides transmitidas . . .
// ------------------------------------------------------------------------------------
string filename = "C:\\Users\\berna\\source\\repos\\TP2\\ifag0780.01n";

// para hacerla fácil, anoto el tiempo correspondiente a las
//  efemérides transmitidas que voy a usar para cada satélite ...
int linea = 0;
string sat = "";
string sata = "";
DateTime fechaHora = new DateTime(2001, 3, 19, 0, 15, 0);
double a0 = 0.0, a1 = 0.0, a2 = 0.0;
double Toe = 0.0, sqrtA = 0.0, Delta_n = 0.0, GPS_Week = 0.0;
Int32 L_sec, GPSweek = 0;
double e = 0.0, M0 = 0.0, i0 = 0.0, idot = 0.0;
double omega = 0.0, OMEGA = 0.0, OMEGA_DOT = 0.0;
double Cuc = 0.0, Cus = 0.0, Crc = 0.0, Crs = 0.0, Cic = 0.0, Cis = 0.0;
double IODE, c2, L2_P, SVa, SVh, TGD, IODC, TxToM, FitInt;
bool start = false;

foreach (string line in System.IO.File.ReadLines(filename))
{
    if (start)
    {
        if (line[5] != '.')
        {
            if (sat != "")
            {
                if ( !satlist.Contains(sat))
                {
                    satlist.Append(sat);
                    mensajes[sata] = new List<Efemerides>();
                    tefem[sata] = fechaHora;
                }
                Efemerides buff = new Efemerides {
                    FechaHora= fechaHora, a0= a0 ,a1= a1,a2= a2,
                    T0e= Toe,GPSweek= GPS_Week,sqrtA= sqrtA,e= e, M0= M0,
                    omega= omega,i0= i0,OMEGA= OMEGA,Delta_n= Delta_n,
                    idot= idot,OMEGA_DOT= OMEGA_DOT,Cuc= Cuc,Cus= Cus,
                    Crc= Crc,Crs= Crs,Cic= Cic,Cis= Cis};
                mensajes[sata].Add(buff);
                GPSweek = Convert.ToInt32(GPS_Week);
            }                        
            sat = line[0..2];
            if (sat.StartsWith(" "))
                sata = "0" + sat[1];
            else
                sata = sat;
            string seg = line[18..22].Replace(" ","0");
            string AA = line[3..5];
            if (Int32.Parse(AA) < 10)
                AA = "0" + line[4];
            //                  AA                         MM                                  DD             
            string Sfecha = "20" + AA +"/"+ line[6..8].Replace(" ", "0") + "/"+ line[9..11].Replace(" ", "0")
            //                           HH                                     mm                    ss
                    + " " + line[12..14].Replace(" ", "0") + ":" + line[15..17].Replace(" ", "0") + ":" + seg;

            CultureInfo provider = CultureInfo.InvariantCulture;
            //DateTime dt = DateTime.ParseExact("2009-05-08 14:40:52,531", "yyyy-MM-dd HH:mm:ss,fff", provider);

            fechaHora = DateTime.ParseExact(Sfecha,"yyyy/MM/dd HH:mm:ss.f",provider);
            a0 = Double.Parse(line[22..41].Replace('D', 'E'));
            a1 = Double.Parse(line[41..60].Replace('D', 'E'));
            a2 = Double.Parse(line[60..79].Replace('D', 'E'));
            linea = 0;
        }
        else if (line.StartsWith("   "))
        {
            linea ++;
            if (linea== 1)
            {
                IODE = Double.Parse(line[3..22].Replace('D', 'E'));
                Crs = Double.Parse(line[22..41].Replace('D', 'E'));
                Delta_n = Double.Parse(line[41..60].Replace('D', 'E'));
                M0 = Double.Parse(line[60..79].Replace('D', 'E'));
            }
            else if (linea== 2)
            {
                Cuc = Double.Parse(line[3..22].Replace('D', 'E'));
                e = Double.Parse(line[22..41].Replace('D', 'E'));
                Cus = Double.Parse(line[41..60].Replace('D', 'E'));
                sqrtA = Double.Parse(line[60..79].Replace('D', 'E'));
            }
            else if (linea== 3)
            {
                Toe = Double.Parse(line[3..22].Replace('D', 'E'));
                Cic = Double.Parse(line[22..41].Replace('D', 'E'));
                OMEGA = Double.Parse(line[41..60].Replace('D', 'E'));
                Cis = Double.Parse(line[60..79].Replace('D', 'E'));
            }
            else if (linea== 4)
            {
                i0 = Double.Parse(line[3..22].Replace('D', 'E'));
                Crc = Double.Parse(line[22..41].Replace('D', 'E'));
                omega = Double.Parse(line[41..60].Replace('D', 'E'));
                OMEGA_DOT = Double.Parse(line[60..79].Replace('D', 'E'));
            }
            else if (linea== 5)
            {                        
                idot = Double.Parse(line[3..22].Replace('D', 'E'));
                c2 = Double.Parse(line[22..41].Replace('D', 'E'));
                GPS_Week = Double.Parse(line[41..60].Replace('D', 'E'));
                L2_P = Double.Parse(line[60..79].Replace('D', 'E'));
            }
            else if (linea== 6)
            {
                SVa = Double.Parse(line[3..22].Replace('D', 'E'));
                SVh = Double.Parse(line[22..41].Replace('D', 'E'));
                TGD = Double.Parse(line[41..60].Replace('D', 'E'));
                IODC = Double.Parse(line[60..79].Replace('D', 'E'));
            }
            else if (linea== 7)
            {
                TxToM = Double.Parse(line[3..22].Replace('D', 'E'));
                FitInt = Double.Parse(line[22..41].Replace('D', 'E'));
            }
        }        
    }
    if (line.Contains("LEAP SECONDS"))
        L_sec = Int32.Parse(line[0..8]);
    if (line.Contains("END OF HEADER"))
        start = true;
}
TimeSpan ts = new TimeSpan(GPSweek * 7, 0, 0, 0);
tGPS0.Add(ts);
// tGPS0 += timedelta(days = GPSweek * 7);
// ------------------------------------------------------------------------------------
//  / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / -
// ------------------------------------------------------------------------------------


// # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


////////////////////////////////////////////////////////////// -  -  -  -  -  -  -  -  -  -  -  -  -
//   I N I C I O    -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
////////////////////////////////////////////////////////////// -  -  -  -  -  -  -  -  -  -  -  -  -

//Console.Clear();

// consigo las coordenadas de cada satélite
//    en el tiempo de emisión de la señal
//    para usarlas en lugar de las precisas
//Calculadas = coord_obs();
/**/
foreach (var s in Precisas.Keys)
{
    int j = 0;
    double[] fila = new double[4];
    foreach (double valor in Precisas[s])
    {
        if (j < 3)
            fila[j] = valor * 1000.0;
        else
            fila[j] = valor;
        j++;
    }
    Calculadas.Add(s, fila);
}
/**/
/*
Console.WriteLine("Efemerides:");
foreach(string st in Calculadas.Keys)
{
    string reng =st + ": ";
    double[] lin = Calculadas[st];
    for (int j=0;j<lin.Length;j++)
            reng += string.Format("{0,20:F16}  ",lin[j]);
    Console.WriteLine(reng);

}
*/


// Muestra las diferencias entre las coordendas precisas del satélite
// al instante de observación y las calculadas al instante de emisión
imprime_dif_sat();

Console.WriteLine(nl+nl+"--------------------------------------------------------");
imprime_Correg();   // Primero muestra la condición inicial desde donde partimos

double Corr_reloj = 0.0;

DateTime inicio = DateTime.Now;
for (int paso=0; paso<1; paso++  )
{
    if (paso == 1)
        inicio = DateTime.Now;

    Console.WriteLine("--------------------------------------------------------");
    Console.WriteLine(string.Format("----> Paso: {0,4:D}",paso + 1));
    Console.WriteLine();

    rr = Coord[0..3];    // Coordenadas (x, Y, Z) calculadas de la estación

    L = np.empty((0, 1));
    A = np.empty((0, 2));
    arma_matriz(ref A, ref L);


    //P = linalg.inv(C)

    //X1 = np.matmul(np.linalg.inv(np.matmul(np.transpose(A), A)), np.matmul(np.transpose(A), L));

    float rcond = float.NaN;
    var XX = np.linalg.lstsq(A, L, rcond);
    Console.WriteLine(nl+"-->> Resultados de la Regresión lineal:");
    Console.WriteLine(XX.Item1.repr);
    Console.WriteLine(XX.Item2.repr);
    Console.WriteLine(XX.Item3);
    Console.WriteLine(XX.Item4.repr);
    X1 = XX.Item1;

    imprime_resu();

    if (paso < 1)
    {
        for (int i = 0; i < 4; i++) Coord[i] = Coord[i] + (double) X1[i];
        Corr_reloj = Coord[3];
    }
    else
    {
        for (int i = 0; i < 3; i++) Coord[i] = Coord[i] + (double) X1[i];
        Coord[3]=Corr_reloj;
    }
    Console.WriteLine(string.Format("Angulo: {0,8:F3}º",calcula_angulo_dif())+nl+nl);
    imprime_Correg();

    //Cxyz = linalg.inv(transpose(A) @ P @ A)
    //Console.WriteLine(Cxyz)
}
DateTime final = DateTime.Now;
double tproceso = (final - inicio).TotalSeconds;
Console.WriteLine(nl+"------------------------------------------------------");
Console.WriteLine(nl+string.Format("tiempo de proceso: {0,8:F5}seg",tproceso)+nl+nl);

/*
Console.Write($"{Environment.NewLine}Press any key to exit...");
Console.ReadKey(true);
*/
public class Efemerides
{
    public DateTime FechaHora { get; set; }
    public double a0 { get; set; }
    public double a1 { get; set; }
    public double a2 { get; set; }
    public double T0e { get; set; }
    public double GPSweek { get; set; }
    public double sqrtA { get; set; }
    public double e { get; set; }
    public double M0 { get; set; }
    public double omega { get; set; }
    public double i0 { get; set; }
    public double OMEGA { get; set; }
    public double Delta_n { get; set; }
    public double idot { get; set; }
    public double OMEGA_DOT { get; set; }
    public double Cuc { get; set; }
    public double Cus { get; set; }
    public double Crc { get; set; }
    public double Crs { get; set; }
    public double Cic { get; set; }
    public double Cis { get; set; }
}
