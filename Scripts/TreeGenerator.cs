using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEditor;



public class TreeGenerator : EditorWindow
{
    public GameObject tree;
    public Material trunkMaterial;
    public Material leavesMaterial;
    public int _recursionLevel = 3;
    public float _roughness = 1f;
    public int _basePolygon = 20;
    public int _floorNumber = 40;
    public float _trunkThickness = 0.35f;
    public float _floorHeight = 0.1f;
    public float _firstBranchHeight = 0.2f;
    public float _distorsionCone = 20f;
    public float _twistiness = 8f;
    public float _reductionRate = 0.1f;
    public float _branchDensity = 0.2f;
    public Mesh _leavesMesh;
    public float _leavesSize = 30f;
    [Range(0, 2)]
    public int smoothingPower = 1;

    private float _randomRatio = 0.005f;
    private Vector2 scrollPosition = new Vector2();

    [MenuItem("Window/TreeGenerator")]
    public static void ShowWindow()
    {

        EditorWindow win = EditorWindow.GetWindow(typeof(TreeGenerator));
        win.titleContent = new GUIContent("TreeGenerator");
        win.minSize = new Vector2(275, 360);
    }

    private void OnGUI()
    {
        if(!trunkMaterial) trunkMaterial = AssetDatabase.LoadAssetAtPath<Material>("Assets/TreeProceduralGeneration/Materials/Trunk.mat");
        if(!leavesMaterial) leavesMaterial = AssetDatabase.LoadAssetAtPath<Material>("Assets/TreeProceduralGeneration/Materials/Leaves.mat");
        if(!tree) tree = GameObject.Find("tree");
        //if(!_leavesMesh) _leavesMesh = Resources.GetBuiltinResource<Mesh>("Sphere.fbx");

        Editor e = Editor.CreateEditor(this);
        scrollPosition = EditorGUILayout.BeginScrollView(scrollPosition, GUIStyle.none, GUI.skin.verticalScrollbar);
        e.DrawDefaultInspector();
        EditorGUILayout.EndScrollView();
        GUILayout.FlexibleSpace();
        if (GUILayout.Button("Generate Tree", GUILayout.Height(64)))
            gen(); 
        if (GUILayout.Button("Save mesh to Assets", GUILayout.Height(64)))
            saveMesh();
    }

    private void saveMesh()
    {
        AssetDatabase.CreateAsset(tree.GetComponent<MeshFilter>().sharedMesh, "Assets/tree.asset");
        AssetDatabase.SaveAssets();
    }

    static T[] SubArray<T>(T[] data, int index, int length)
    {
        T[] result = new T[length];
        for (int i = index; i < index + length; ++i)
        {
            result[i] = data[i];
        }
        return result;
    }

    int getCloseValue(int val)
    {
        return Random.Range(val - (int)(val * _randomRatio), val + (int)(val * _randomRatio));
    }

    float getCloseValue(float val)
    {
        return Random.Range(val - val * _randomRatio, val + val * _randomRatio);
    }

    // Start is called before the first frame update
    void gen()
    {
        if (!_leavesMesh)
            _leavesMesh = new Mesh();
        if (!leavesMaterial)
            leavesMaterial = trunkMaterial;
        if (tree)
            DestroyImmediate(tree);
        tree = new GameObject("tree");
        

        int basePolygon = Mathf.Max(3, getCloseValue(_basePolygon));
        Vector3[] startVertices = new Vector3[basePolygon];
        Vector2[] startUv = new Vector2[basePolygon];

        float angularStep = 2f * Mathf.PI / (float)(basePolygon);
        for (int j = 0; j < basePolygon; ++j)
        {
            Vector3 randomness = new Vector3
                (
                Random.Range(-_roughness, _roughness),
                Random.Range(-_roughness, _roughness),
                Random.Range(-_roughness, _roughness)
                ) / 10f;
            Vector3 pos = new Vector3(Mathf.Cos(j * angularStep), 0f, Mathf.Sin(j * angularStep)) + randomness;
            startVertices[j] = _trunkThickness * (pos);
            startUv[j] = new Vector2(j * angularStep, startVertices[j].y);
        }

        Mesh mesh = GenBranch(tree, basePolygon, startVertices, startUv, getCloseValue(_trunkThickness), getCloseValue(_floorHeight), getCloseValue(_floorNumber), new Vector3(), Vector3.up, 0f, getCloseValue(_distorsionCone), getCloseValue(_roughness), getCloseValue(_branchDensity), getCloseValue(_recursionLevel));

        CombineInstance[] combine = new CombineInstance[2];
        combine[0].mesh = MeshSmoothener.SmoothMesh(mesh.GetSubmesh(0), smoothingPower, MeshSmoothener.Filter.Laplacian);
        combine[0].transform = tree.transform.localToWorldMatrix;
        combine[1].mesh = mesh.GetSubmesh(1);
        combine[1].transform = tree.transform.localToWorldMatrix;
        mesh.CombineMeshes(combine, false, false);
        //mesh = MeshSmoothener.SmoothMesh(mesh, smoothingPower, MeshSmoothener.Filter.Laplacian);

        mesh.RecalculateNormals();

        tree.AddComponent<MeshFilter>().mesh = mesh;
        tree.AddComponent<MeshRenderer>().materials = new Material[2] { trunkMaterial, leavesMaterial };

    }

    Vector3 ChangeCoordinates(Vector3 input, Vector3 inputNormal, Vector3 newNormal)
    {
        float angle = Vector3.Angle(inputNormal, newNormal);
        Vector3 axis = Vector3.Cross(inputNormal, newNormal);
        Quaternion rot = Quaternion.AngleAxis(angle, axis);
        return rot * input;
    }

    Vector3 getRandomVectorInCone(float coneAngularAmplitude, Vector3 direction)
    {
        return (new Vector3(Random.Range(-coneAngularAmplitude, coneAngularAmplitude), Random.Range(-coneAngularAmplitude, coneAngularAmplitude), Random.Range(-coneAngularAmplitude, coneAngularAmplitude)) / 100f + direction).normalized;
    }

    static bool Happens(float proba)
    {
        return Random.Range(0f, 1f) < proba;
    }
    
    Mesh GenBranch(GameObject tree, int basePolygon, Vector3[] startVertices, Vector2[] startUv, float thickness, float floorHeight, int floorNumber, Vector3 startingPos, Vector3 startingDirection, float angularOffset, float distorsionCone, float roughness, float branchDensity, int recursionLevel)
    {
        float currentThickness = thickness;
        Mesh mesh = new Mesh();
        mesh.subMeshCount = 2;

        Vector3[] vertices = new Vector3[basePolygon * floorNumber];
        Vector2[] uv = new Vector2[vertices.Length];
        int[] triangles;

        Vector3 leavesPosition = new Vector3();

        float angularStep = 2f * Mathf.PI / (float)(basePolygon);
        Vector3 first = new Vector3(Mathf.Cos(angularOffset), 0f, Mathf.Sin(angularOffset)) * thickness;
        first = ChangeCoordinates(first, Vector3.up, startingDirection);
        first += startingPos;
        
        for (int i = 0; i < startVertices.Length; ++i)
        {
            vertices[i] = startVertices[i];
            uv[i] = startUv[i];
        }
            

        triangles = new int[6 * (vertices.Length - basePolygon)];


        Vector3 growDirection = startingDirection;

        Vector3 lastPivot = startingPos;

        for (int i = 1; i < vertices.Length/basePolygon; ++i)
        {
            
            Vector3 pivot = lastPivot + floorHeight * growDirection;
            lastPivot = pivot;

            for (int j = 0; j < basePolygon; ++j)
            {
                Vector3 randomness = new Vector3
                    (
                    Random.Range(-roughness / 10f, roughness / 10f),
                    Random.Range(-roughness / 10f, roughness / 10f),
                    Random.Range(-roughness / 10f, roughness / 10f)
                    );
                Vector3 pos = new Vector3(Mathf.Cos(j * angularStep + angularOffset), 0f, Mathf.Sin(j * angularStep + angularOffset)) + randomness;
                pos *= currentThickness;
                pos = ChangeCoordinates(pos, new Vector3(0f, 1f, 0f), growDirection);
                vertices[i * basePolygon + j] = pos + pivot;
                uv[i * basePolygon + j] = new Vector2( j*angularStep, vertices[i * basePolygon + j].y );
                
                triangles[6 * ((i - 1) * basePolygon + j)]     = (i - 1) * basePolygon + j;
                triangles[6 * ((i - 1) * basePolygon + j) + 1] = (i) * basePolygon + j;
                triangles[6 * ((i - 1) * basePolygon + j) + 2] = (i - 1) * basePolygon + (j + 1) % basePolygon;
                triangles[6 * ((i - 1) * basePolygon + j) + 3] = (i - 1) * basePolygon + (j + 1) % basePolygon;
                triangles[6 * ((i - 1) * basePolygon + j) + 4] = (i) * basePolygon + j;
                triangles[6 * ((i - 1) * basePolygon + j) + 5] = (i) * basePolygon + (j + 1) % basePolygon;
            }
            if (Happens(branchDensity) && recursionLevel > 0 && i >= floorNumber * _firstBranchHeight && basePolygon >= 4 && i < vertices.Length / basePolygon - 1) // split!
            {
                int subBasePolygon = Random.Range(Mathf.Max(2, basePolygon/3), Mathf.Min(2*basePolygon/3, basePolygon - 1)) + 1;
                int subBasePolygon2 = basePolygon - subBasePolygon + 2;
                int newOffset = Random.Range(0, basePolygon);
                Vector3[] subStartVertices = new Vector3[subBasePolygon];
                Vector2[] subStartUv = new Vector2[subBasePolygon];
                Vector3[] subStartVertices2 = new Vector3[subBasePolygon2];
                Vector2[] subStartUv2 = new Vector2[subBasePolygon2];

                Vector3 subStartingPos = pivot;
                Vector3 mid = new Vector3();
                for (int k = 0; k < subBasePolygon; ++k)
                {
                    int shift = ((k + newOffset) % basePolygon + basePolygon) % basePolygon;
                    subStartVertices[k] = vertices[i * basePolygon + shift];
                    subStartUv[k] = uv[i * basePolygon + shift];
                    mid += subStartVertices[k];
                }

                mid /= subBasePolygon;
                subStartingPos = mid;

                float newAngularOffset = angularOffset + Vector3.SignedAngle(subStartVertices[0] - mid, vertices[i * basePolygon] - pivot, growDirection) * Mathf.Deg2Rad;


                Vector3 subStartingPos2 = pivot;
                Vector3 mid2 = new Vector3();
                for (int k = 0; k < subBasePolygon2; ++k)
                {
                    int shift = ((k + newOffset + subBasePolygon - 1) % basePolygon + basePolygon) % basePolygon;
                    subStartVertices2[k] = vertices[i * basePolygon + shift];
                    subStartUv2[k] = uv[i * basePolygon + shift];
                    mid2 += subStartVertices2[k];
                }
                mid2 /= subBasePolygon2;
                subStartingPos2 = mid2;

                float newAngularOffset2 = angularOffset + Vector3.SignedAngle(subStartVertices2[0] - mid2, vertices[i * basePolygon] - pivot, growDirection) * Mathf.Deg2Rad;

                mesh.vertices = SubArray(vertices, 0, basePolygon * (i + 1));
                mesh.uv = SubArray(uv, 0, basePolygon * (i + 1));
                mesh.SetTriangles(SubArray(triangles, 0, 6 * (i * basePolygon)), 0);
                mesh.RecalculateNormals();

                Mesh res1 = new Mesh();
                res1.subMeshCount = 2;
                Mesh res2 = new Mesh();
                res2.subMeshCount = 2;
                Mesh res = new Mesh();
                res.subMeshCount = 2;

                Mesh branch1 = GenBranch(tree, subBasePolygon,  subStartVertices,  subStartUv, 
                    currentThickness * subBasePolygon  / ((float)basePolygon), floorHeight, 
                    floorNumber - i, subStartingPos,  getRandomVectorInCone(distorsionCone, (growDirection + distorsionCone / 45f * (mid - pivot).normalized).normalized), 
                    newAngularOffset, distorsionCone, roughness, branchDensity * 1.1f, recursionLevel - 1);

                Mesh branch2 = GenBranch(tree, subBasePolygon2, subStartVertices2, subStartUv2,
                    currentThickness * subBasePolygon2 / ((float)basePolygon), floorHeight,
                    floorNumber - i, subStartingPos2, getRandomVectorInCone(distorsionCone, (growDirection + distorsionCone / 45f * (mid2 - pivot).normalized).normalized),
                    newAngularOffset2, distorsionCone, roughness, branchDensity * 1.1f, recursionLevel - 1);
                
                CombineInstance[] combine = new CombineInstance[3];
                combine[1].transform = tree.transform.localToWorldMatrix;
                combine[2].transform = tree.transform.localToWorldMatrix;

                Mesh subMesh1 = mesh.GetSubmesh(0);
                Mesh subMesh11 = branch1.GetSubmesh(0);
                Mesh subMesh12 = branch2.GetSubmesh(0);

                Mesh subMesh2 = mesh.GetSubmesh(1);
                if (subMesh2 == null)
                {
                    subMesh2 = new Mesh();
                }

                Mesh subMesh21 = branch1.GetSubmesh(1);
                if(subMesh21 == null)
                {
                    subMesh21 = new Mesh();
                }

                Mesh subMesh22 = branch2.GetSubmesh(1);
                if (subMesh22 == null)
                {
                    subMesh22 = new Mesh();
                }

                combine[0].mesh = subMesh1;
                combine[0].transform = tree.transform.localToWorldMatrix;
                combine[1].mesh = subMesh11;
                combine[1].transform = tree.transform.localToWorldMatrix;
                combine[2].mesh = subMesh12;
                combine[2].transform = tree.transform.localToWorldMatrix;

                res1.CombineMeshes(combine, true, false);
                res1.RecalculateNormals();

                combine[0].mesh = subMesh2;
                combine[1].mesh = subMesh21;
                combine[2].mesh = subMesh22;

                res2.CombineMeshes(combine, true, false);
                res2.RecalculateNormals();

                combine = new CombineInstance[2];
                combine[0].mesh = res1;
                combine[0].transform = tree.transform.localToWorldMatrix;
                combine[1].mesh = res2;
                combine[1].transform = tree.transform.localToWorldMatrix;

                res.CombineMeshes(combine, false, false);
                res.RecalculateNormals();
                res.Optimize();
                return res;
            }
            currentThickness = Mathf.Pow(Mathf.Asin((floorNumber - i) / (float)(floorNumber) * 2f - 1f) / Mathf.PI + 0.5f, _reductionRate) * thickness;
            growDirection = getRandomVectorInCone(_twistiness, growDirection);
            
            leavesPosition = pivot;

        }

        mesh.vertices = vertices;
        mesh.uv = uv;
        mesh.SetTriangles(triangles, 0);
        mesh.RecalculateNormals();
        {

            Mesh trLeavesMesh = new Mesh();
            trLeavesMesh.subMeshCount = 2;

            Vector3[] trVertices = (Vector3[])(_leavesMesh.vertices.Clone());
            int[] leavesTriangles = (int[])(_leavesMesh.triangles.Clone());

            for (int i = 0; i < _leavesMesh.vertices.Length; ++i)
            {
                trVertices[i] *= currentThickness * _leavesSize;
                trVertices[i] += leavesPosition;
            }
            trLeavesMesh.SetVertices(trVertices);
            trLeavesMesh.SetTriangles(leavesTriangles, 0);
            trLeavesMesh.SetNormals((Vector3[])(_leavesMesh.normals.Clone()));
            trLeavesMesh.SetTangents((Vector4[])(_leavesMesh.tangents.Clone()));
            trLeavesMesh.uv = (Vector2[])(_leavesMesh.uv.Clone());

            Mesh res = new Mesh();
            res.subMeshCount = 2;

            CombineInstance[] combine = new CombineInstance[2];
            combine[0].mesh = mesh;
            combine[0].transform = tree.transform.localToWorldMatrix;
            combine[1].mesh = trLeavesMesh;
            combine[1].transform = tree.transform.localToWorldMatrix;
            res.CombineMeshes(combine, false, false);
            res.RecalculateNormals();
            res.Optimize();
            res.OptimizeIndexBuffers();
            res.OptimizeReorderVertexBuffer();
            return res;
        }
    }
}
