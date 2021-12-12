using System;
using System.Collections.Generic;
using System.Numerics;
using UnityEngine;
using UnityEngine.UI;
using Matrix4x4 = UnityEngine.Matrix4x4;
using Quaternion = UnityEngine.Quaternion;
using Vector3 = UnityEngine.Vector3;

public class MySphere
{
    public Vector3 Force;
    public double Pressure;
    public int ObjectPosition_1;
    public int ObjectPosition_2;
    public Vector3 Speed;
    public Vector3 location;
    public PhysObject Parent;
    public MySphere()
    {
        Pressure = 0;
    }
}

public class PhysObject
{
    public MySphere obj_;
    public MySphere left;
    public MySphere right;
    public MySphere front;
    public MySphere back;
    public MySphere up;
    public MySphere down;
    public float stepToRight;
    public float pressure;
    public PhysObject(MySphere new_obj)
    {
        obj_ = new_obj;
        obj_.Parent = this;
    }
}
public class Problem : MonoBehaviour
{
    // Start is called before the first frame update
    private List<Transform> InstantiateObjects;
    public int XObjectNumber;
    public float XObjectStep;
    public int YObjectNumber;
    public float YObjectStep;
    public int ZObjectNumber;
    public float ZObjectStep;
    public float K = 10;
    public float Mass = 1;
    public float C = 20;
    public Text PotentialEnergy;
    public double Pe;
    public Text KineticEnergy;
    public double Ke;
    public Text SummEnergy;
    public double Summ;
    private List<Matrix4x4[]> matrices = new List<Matrix4x4[]>();
    private Mesh mesh;
    private int population;
    private MaterialPropertyBlock block;
    public Material material;
    private PhysObject[][][] ObjectsArray;
    private float distance;
    // Plot
    public SimplestPlot.PlotType PlotExample = SimplestPlot.PlotType.TimeSeries;
    public int DataPoints = 10000;
    public SimplestPlot EnergyPlot;
    public SimplestPlot PressurePlot;
    private float Counter = 0;
    private Color[] MyColors = new Color[3];

    private float[] XValues; // just X
    private float[] PEValues; // Potential energy
    private float[] KEValues; // Kinetic energy
    private float[] SEValues; // Full energy
    private float[] Pressure; // Full energy
    private float[] CountSamePres; // Full energy
    private float SpeedOfProg;
    private bool is_plot_not_defined = true;
    void Start()
    {
        population = XObjectNumber * YObjectNumber * ZObjectNumber;
        for (int i = 0; i < population / 1000; ++i)
        {
            var matrices_temp = new Matrix4x4[1000];
            matrices.Add(matrices_temp);
        }
        matrices.Add(new Matrix4x4[population % 1000]);
        var Object = GameObject.CreatePrimitive(PrimitiveType.Sphere);
        Object.active = false;
        mesh = Object.GetComponent<MeshFilter>().mesh;
        block = new MaterialPropertyBlock();
        //create objects
        ObjectCreating();
        //connect objects
        ObjectConnecting();
        //
        SetForce();
        distance = (matrices[0][0].GetColumn(3) - matrices[0][matrices[0].Length - 1].GetColumn(3)).magnitude;
        var deltaX = (matrices[0][0].GetColumn(3) - matrices[0][matrices[0].Length - 1].GetColumn(3)).magnitude;
        var RestoringForce = -(distance - deltaX) * K;
        PotentialEnergy.text = "Potential energy = " + Pe;
        KineticEnergy.text = "Kinetic energy = " + Ke;
        SummEnergy.text = "Summ = " + (Ke + Pe);
        //EnergyPlot
        XValues = new float[DataPoints];
        PEValues = new float[DataPoints];
        KEValues = new float[DataPoints];
        SEValues = new float[DataPoints];
        for (int i = 0; i < DataPoints; ++i)
        {
            XValues[i] = i;
        }
        PEValues[DataPoints - 1] = (float)Pe;
        KEValues[DataPoints - 1] = (float)Ke;
        SEValues[DataPoints - 1] = (float)(Ke + Pe);
        MyColors[0] = Color.red; // Kinetic
        MyColors[1] = Color.green; // Potential
        MyColors[2] = Color.blue; // Summ
        EnergyPlot.SetResolution(new UnityEngine.Vector2(500, 500));
        EnergyPlot.BackGroundColor = new Color(0.1f, 0.1f, 0.1f, 0.4f);
        EnergyPlot.TextColor = Color.yellow;
        // PressurePlot
        Pressure = new float[(XObjectNumber - 1) * YObjectNumber * ZObjectNumber];
        CountSamePres = new float[(XObjectNumber - 1) * YObjectNumber * ZObjectNumber];
        for (int i = 0; i < (XObjectNumber - 1) * YObjectNumber * ZObjectNumber; ++i)
        {
            CountSamePres[i] = i;
        }
        PressurePlot.SetResolution(new UnityEngine.Vector2(500, 500));
        PressurePlot.BackGroundColor = new Color(0.1f, 0.1f, 0.1f, 0.9f);
        PressurePlot.TextColor = Color.yellow;
    }

    void Update()
    {
        if (is_plot_not_defined)
        {
            for (int Cnt = 0; Cnt < 3; ++Cnt)
            {
                EnergyPlot.SeriesPlotY.Add(new SimplestPlot.SeriesClass());
                EnergyPlot.DistributionPlot.Add(new SimplestPlot.DistributionClass());
                EnergyPlot.PhaseSpacePlot.Add(new SimplestPlot.PhaseSpaceClass());
                EnergyPlot.SeriesPlotY[Cnt].MyColor = MyColors[Cnt];
                EnergyPlot.DistributionPlot[Cnt].MyColor = MyColors[Cnt];
                EnergyPlot.PhaseSpacePlot[Cnt].MyColor = MyColors[Cnt];
                EnergyPlot.DistributionPlot[Cnt].NumberOfBins = (Cnt + 1) * 5;
            }
            PressurePlot.SeriesPlotY.Add(new SimplestPlot.SeriesClass());
            PressurePlot.DistributionPlot.Add(new SimplestPlot.DistributionClass());
            PressurePlot.PhaseSpacePlot.Add(new SimplestPlot.PhaseSpaceClass());
            PressurePlot.SeriesPlotY[0].MyColor = MyColors[0];
            PressurePlot.DistributionPlot[0].MyColor = MyColors[0];
            PressurePlot.PhaseSpacePlot[0].MyColor = MyColors[0];
            is_plot_not_defined = false;
        }
        for (int i = 0; i < matrices.Count; ++i)
            Graphics.DrawMeshInstanced(mesh, 0, material, matrices[i], matrices[i].Length, block);
        for (int i = 0; i < XObjectNumber - 1; ++i)
        {
            for (int j = 0; j < YObjectNumber; ++j)
            {
                for (int k = 0; k < ZObjectNumber; ++k)
                {
                    if (ObjectsArray[i][j][k] != null && ObjectsArray[i + 1][j][k] != null)
                        Debug.DrawLine(ObjectsArray[i][j][k].obj_.location, ObjectsArray[i + 1][j][k].obj_.location);
                    //if (j + 1 < YObjectNumber && k + 1 < ZObjectNumber && ObjectsArray[i][j][k] != null && ObjectsArray[i][j + 1][k + 1] != null)
                    //Debug.DrawLine(ObjectsArray[i][j][k].obj_.location, ObjectsArray[i][j + 1][k + 1].obj_.location);
                    if (j + 1 < YObjectNumber && ObjectsArray[i][j][k] != null && ObjectsArray[i][j + 1][k] != null)
                        Debug.DrawLine(ObjectsArray[i][j][k].obj_.location, ObjectsArray[i][j + 1][k].obj_.location);
                    if (k + 1 < ZObjectNumber && ObjectsArray[i][j][k] != null && ObjectsArray[i][j][k + 1] != null)
                        Debug.DrawLine(ObjectsArray[i][j][k].obj_.location, ObjectsArray[i][j][k + 1].obj_.location);

                }
            }
        }
        //
        Counter++;
        EnergyPlot.MyPlotType = PlotExample;
        switch (PlotExample)
        {
            case SimplestPlot.PlotType.TimeSeries:
                EnergyPlot.SeriesPlotY[0].YValues = PEValues;
                EnergyPlot.SeriesPlotY[1].YValues = KEValues;
                EnergyPlot.SeriesPlotY[2].YValues = SEValues;
                EnergyPlot.SeriesPlotX = XValues;
                PressurePlot.SeriesPlotY[0].YValues = Pressure;
                PressurePlot.SeriesPlotX = CountSamePres;
                break;
            case SimplestPlot.PlotType.Distribution:
                EnergyPlot.DistributionPlot[0].Values = PEValues;
                EnergyPlot.DistributionPlot[1].Values = KEValues;
                EnergyPlot.DistributionPlot[2].Values = SEValues;
                break;
            case SimplestPlot.PlotType.PhaseSpace:
                EnergyPlot.PhaseSpacePlot[0].XValues = XValues;
                EnergyPlot.PhaseSpacePlot[0].YValues = PEValues;
                EnergyPlot.PhaseSpacePlot[1].XValues = PEValues;
                EnergyPlot.PhaseSpacePlot[1].YValues = KEValues;

                break;
            default:
                break;
        }
        EnergyPlot.UpdatePlot();
        PressurePlot.UpdatePlot();
        UpdateS();
    }
    // Update is called once per frame
    //void FixedUpdate()
    void UpdateS()
    {
        long summTime = 0;
        var watch = System.Diagnostics.Stopwatch.StartNew();
        while (summTime < 10)
        {
            var delta = 0.0001f;
            for (int i = 0; i < XObjectNumber - 1; ++i)
            {
                for (int j = 0; j < YObjectNumber; ++j)
                {
                    for (int k = 0; k < ZObjectNumber; ++k)
                    {
                        var obj = ObjectsArray[i][j][k];
                        if (obj == null) continue;
                        obj.obj_.Speed += obj.obj_.Force * (delta / Mass);
                        Ke += Math.Pow(obj.obj_.Speed.magnitude, 2) * Mass * 0.5;
                        var t = obj.obj_.Speed * delta; // x = v * t
                        obj.obj_.location += t;
                        var m = Matrix4x4.Translate(t);
                        matrices[obj.obj_.ObjectPosition_2][obj.obj_.ObjectPosition_1] =
                            Matrix4x4.TRS(obj.obj_.location, Quaternion.identity, new Vector3(1, 1, 1));
                    }
                }
            }

            for (int i = 0; i < XObjectNumber - 1; ++i)
            {
                for (int j = 0; j < YObjectNumber; ++j)
                {
                    for (int k = 0; k < ZObjectNumber; ++k)
                    {
                        var obj = ObjectsArray[i][j][k];
                        if (obj == null) continue;
                        obj.obj_.Force = Vector3.zero;
                        if (obj.right != null)
                        {
                            var diff = obj.right.location - obj.obj_.location;
                            float x = diff.magnitude - obj.stepToRight;
                            obj.obj_.Force += K * x * diff.normalized;
                            Pe += Math.Pow(x, 2) * 0.5 * K;
                            if (obj.right.Force.magnitude != 0.0f)
                            {
                                obj.pressure = (obj.right.Force - obj.obj_.Force).magnitude;
                            }
                        }

                        if (obj.left != null)
                        {
                            var diff = obj.left.location - obj.obj_.location;
                            float x = diff.magnitude - XObjectStep;
                            obj.obj_.Force += K * x * diff.normalized;
                            //Pe += Math.Pow(x, 2) * 0.5 * K;
                        }

                        if (obj.front != null)
                        {
                            var diff = obj.front.location - obj.obj_.location;
                            float y = diff.magnitude - YObjectStep;
                            obj.obj_.Force += K * y * diff.normalized;
                            Pe += Math.Pow(y, 2) * 0.5 * K;
                        }

                        if (obj.back != null)
                        {
                            var diff = obj.back.location - obj.obj_.location;
                            float y = diff.magnitude - YObjectStep;
                            obj.obj_.Force += K * y * diff.normalized;
                            //Pe += Math.Pow(y, 2) * 0.5 * K;
                        }

                        if (obj.down != null)
                        {
                            var diff = obj.down.location - obj.obj_.location;
                            float z = diff.magnitude - ZObjectStep;
                            obj.obj_.Force += K * z * diff.normalized;
                            //Pe += Math.Pow(z, 2) * 0.5 * K;
                        }

                        if (obj.up != null)
                        {
                            var diff = obj.up.location - obj.obj_.location;
                            float z = diff.magnitude - ZObjectStep;
                            obj.obj_.Force += K * z * diff.normalized;
                            Pe += Math.Pow(z, 2) * 0.5 * K;
                        }
                    }
                }
            }

            watch.Stop();
            summTime += watch.ElapsedMilliseconds;
            watch.Start();
        }

        Pe *= 200;
        Ke *= 200;
        PotentialEnergy.text = "Pot/energy = " + Pe;
        KineticEnergy.text = "Kin/energy = " + Ke;
        SummEnergy.text = "Summ = " + (Ke + Pe);
        for (int i = 0; i < DataPoints - 1; ++i)
        {
            PEValues[i] = PEValues[i + 1];
            KEValues[i] = KEValues[i + 1];
            SEValues[i] = SEValues[i + 1];
        }

        PEValues[DataPoints - 1] = (float)Pe;
        KEValues[DataPoints - 1] = (float)Ke;
        SEValues[DataPoints - 1] = (float)(Ke + Pe);
        Pe = 0;
        Ke = 0;
        int counter = 0;
        for (int i = 0; i < XObjectNumber - 1; ++i)
        {
            for (int j = 0; j < YObjectNumber; ++j)
            {
                for (int k = 0; k < ZObjectNumber; ++k)
                {
                    if (ObjectsArray[i][j][k] != null)
                        Pressure[counter++] = ObjectsArray[i][j][k].pressure;
                }
            }
        }
    }

    void SetForce()
    {
        //for (int i = 0; i < XObjectNumber - 1; ++i)
        for (int j = 0; j < YObjectNumber; ++j)
            for (int k = 0; k < ZObjectNumber; ++k)
            {
                if (ObjectsArray[0][j][k] != null)
                    ObjectsArray[0][j][k].obj_.Force = Vector3.right * 5000;
                //ObjectsArray[1][0][0].obj_.Force = Vector3.right * 100;
                //ObjectsArray[0][0][1].obj_.Force = Vector3.right * 1;
            }
    }

    void ObjectCreating()
    {
        ObjectsArray = new PhysObject[XObjectNumber][][];
        for (int i = 0; i < ObjectsArray.Length; ++i)
        {
            ObjectsArray[i] = new PhysObject[YObjectNumber][];
            for (int j = 0; j < ObjectsArray[i].Length; ++j)
            {
                ObjectsArray[i][j] = new PhysObject[ZObjectNumber];
            }
        }
        int pos = 0;
        int pos_2 = 0;
        for (int i = 0; i < XObjectNumber; ++i)
        {
            for (int j = 0; j < YObjectNumber; ++j)
            {
                for (int k = 0; k < ZObjectNumber; ++k)
                {
                    //for empty center
                    //if (i > 0 && k > 0 && j > 0 && i < XObjectNumber - 1 && j < YObjectNumber - 1 && k < ZObjectNumber - 1) continue;
                    ObjectsArray[i][j][k] = new PhysObject(new MySphere());
                    var objMatr = new Matrix4x4();
                    ObjectsArray[i][j][k].obj_.location =
                        new Vector3(i * XObjectStep, j * YObjectStep, k * ZObjectStep);

                    objMatr.SetTRS(ObjectsArray[i][j][k].obj_.location,
                        Quaternion.identity,
                        new Vector3(2, 2, 2));

                    if (pos < 1000)
                    {
                        ObjectsArray[i][j][k].obj_.ObjectPosition_1 = pos;
                        ObjectsArray[i][j][k].obj_.ObjectPosition_2 = pos_2;
                        matrices[pos_2][pos++] = objMatr;
                    }
                    else
                    {
                        pos = 0;
                        ++pos_2;
                    }
                }
            }
        }
    }

    void ObjectConnecting()
    {
        for (int i = 0; i < XObjectNumber; ++i)
        {
            for (int j = 0; j < YObjectNumber; ++j)
            {
                for (int k = 0; k < ZObjectNumber; ++k)
                {
                    if (ObjectsArray[i][j][k] == null) continue;
                    if (i > 0 && ObjectsArray[i - 1][j][k] != null)
                    {
                        ObjectsArray[i][j][k].left = ObjectsArray[i - 1][j][k].obj_;
                    }
                    if (i < XObjectNumber - 1 && ObjectsArray[i + 1][j][k] != null)
                    {
                        ObjectsArray[i][j][k].right = ObjectsArray[i + 1][j][k].obj_;
                        ObjectsArray[i][j][k].stepToRight = (i + 1 - i) * XObjectStep;
                    }
                    if (j > 0 && ObjectsArray[i][j - 1][k] != null)
                    {
                        ObjectsArray[i][j][k].back = ObjectsArray[i][j - 1][k].obj_;
                    }
                    if (j < YObjectNumber - 1 && ObjectsArray[i][j + 1][k] != null)
                    {
                        ObjectsArray[i][j][k].front = ObjectsArray[i][j + 1][k].obj_;
                    }
                    if (k > 0 && ObjectsArray[i][j][k - 1] != null)
                    {
                        ObjectsArray[i][j][k].down = ObjectsArray[i][j][k - 1].obj_;
                    }
                    if (k < ZObjectNumber - 1 && ObjectsArray[i][j][k + 1] != null)
                    {
                        ObjectsArray[i][j][k].up = ObjectsArray[i][j][k + 1].obj_;
                    }

                    if (i != -1 && ObjectsArray[i][j][k].right == null)
                    {
                        ObjectsArray[i][j][k].right = ObjectsArray[XObjectNumber - 1][j][k].obj_;
                        ObjectsArray[i][j][k].stepToRight = (XObjectNumber - 1 - i) * XObjectStep;
                    }
                }
            }
        }
    }
}
