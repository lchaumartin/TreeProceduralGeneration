using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class MeshUtils
{
	// Finds a set of adjacent vertices for a given vertex
	// Note the success of this routine expects only the set of neighboring faces to eacn contain one vertex corresponding
	// to the vertex in question
	public static List<Vector3> findAdjacentNeighbors(Vector3[] v, int[] t, Vector3 vertex)
	{
		List<Vector3> adjacentV = new List<Vector3>();
		List<int> facemarker = new List<int>();
		int facecount = 0;

		// Find matching vertices
		for (int i = 0; i < v.Length; i++)
			if (Mathf.Approximately(vertex.x, v[i].x) &&
				Mathf.Approximately(vertex.y, v[i].y) &&
				Mathf.Approximately(vertex.z, v[i].z))
			{
				int v1 = 0;
				int v2 = 0;
				bool marker = false;

				// Find vertex indices from the triangle array
				for (int k = 0; k < t.Length; k = k + 3)
					if (facemarker.Contains(k) == false)
					{
						v1 = 0;
						v2 = 0;
						marker = false;

						if (i == t[k])
						{
							v1 = t[k + 1];
							v2 = t[k + 2];
							marker = true;
						}

						if (i == t[k + 1])
						{
							v1 = t[k];
							v2 = t[k + 2];
							marker = true;
						}

						if (i == t[k + 2])
						{
							v1 = t[k];
							v2 = t[k + 1];
							marker = true;
						}

						facecount++;
						if (marker)
						{
							// Once face has been used mark it so it does not get used again
							facemarker.Add(k);

							// Add non duplicate vertices to the list
							if (isVertexExist(adjacentV, v[v1]) == false)
							{
								adjacentV.Add(v[v1]);
								//Debug.Log("Adjacent vertex index = " + v1);
							}

							if (isVertexExist(adjacentV, v[v2]) == false)
							{
								adjacentV.Add(v[v2]);
								//Debug.Log("Adjacent vertex index = " + v2);
							}
							marker = false;
						}
					}
			}

		//Debug.Log("Faces Found = " + facecount);

		return adjacentV;
	}


	// Finds a set of adjacent vertices indexes for a given vertex
	// Note the success of this routine expects only the set of neighboring faces to eacn contain one vertex corresponding
	// to the vertex in question
	public static List<int> findAdjacentNeighborIndexes(Vector3[] v, int[] t, Vector3 vertex)
	{
		List<int> adjacentIndexes = new List<int>();
		List<Vector3> adjacentV = new List<Vector3>();
		List<int> facemarker = new List<int>();
		int facecount = 0;

		// Find matching vertices
		for (int i = 0; i < v.Length; i++)
			if (Mathf.Approximately(vertex.x, v[i].x) &&
				Mathf.Approximately(vertex.y, v[i].y) &&
				Mathf.Approximately(vertex.z, v[i].z))
			{
				int v1 = 0;
				int v2 = 0;
				bool marker = false;

				// Find vertex indices from the triangle array
				for (int k = 0; k < t.Length; k = k + 3)
					if (facemarker.Contains(k) == false)
					{
						v1 = 0;
						v2 = 0;
						marker = false;

						if (i == t[k])
						{
							v1 = t[k + 1];
							v2 = t[k + 2];
							marker = true;
						}

						if (i == t[k + 1])
						{
							v1 = t[k];
							v2 = t[k + 2];
							marker = true;
						}

						if (i == t[k + 2])
						{
							v1 = t[k];
							v2 = t[k + 1];
							marker = true;
						}

						facecount++;
						if (marker)
						{
							// Once face has been used mark it so it does not get used again
							facemarker.Add(k);

							// Add non duplicate vertices to the list
							if (isVertexExist(adjacentV, v[v1]) == false)
							{
								adjacentV.Add(v[v1]);
								adjacentIndexes.Add(v1);
								//Debug.Log("Adjacent vertex index = " + v1);
							}

							if (isVertexExist(adjacentV, v[v2]) == false)
							{
								adjacentV.Add(v[v2]);
								adjacentIndexes.Add(v2);
								//Debug.Log("Adjacent vertex index = " + v2);
							}
							marker = false;
						}
					}
			}

		//Debug.Log("Faces Found = " + facecount);

		return adjacentIndexes;
	}

	// Does the vertex v exist in the list of vertices
	static bool isVertexExist(List<Vector3> adjacentV, Vector3 v)
	{
		bool marker = false;
		foreach (Vector3 vec in adjacentV)
			if (Mathf.Approximately(vec.x, v.x) && Mathf.Approximately(vec.y, v.y) && Mathf.Approximately(vec.z, v.z))
			{
				marker = true;
				break;
			}

		return marker;
	}
}

public class SmoothFilter
{
	/*
		Standard Laplacian Smooth Filter
	*/
	public static Vector3[] laplacianFilter(Vector3[] sv, int[] t)
	{
		Vector3[] wv = new Vector3[sv.Length];
		List<Vector3> adjacentVertices = new List<Vector3>();

		float dx = 0.0f;
		float dy = 0.0f;
		float dz = 0.0f;

		for (int vi = 0; vi < sv.Length; vi++)
		{
			// Find the sv neighboring vertices
			adjacentVertices = MeshUtils.findAdjacentNeighbors(sv, t, sv[vi]);

			if (adjacentVertices.Count != 0)
			{
				dx = 0.0f;
				dy = 0.0f;
				dz = 0.0f;

				//Debug.Log("Vertex Index Length = "+vertexIndexes.Length);
				// Add the vertices and divide by the number of vertices
				for (int j = 0; j < adjacentVertices.Count; j++)
				{
					dx += adjacentVertices[j].x;
					dy += adjacentVertices[j].y;
					dz += adjacentVertices[j].z;
				}

				wv[vi].x = dx / adjacentVertices.Count;
				wv[vi].y = dy / adjacentVertices.Count;
				wv[vi].z = dz / adjacentVertices.Count;
			}
		}

		return wv;
	}

	/*
		HC (Humphrey’s Classes) Smooth Algorithm - Reduces Shrinkage of Laplacian Smoother
 
		Where sv - original points
				pv - previous points,
				alpha [0..1] influences previous points pv, e.g. 0
				beta  [0..1] e.g. > 0.5
	*/
	public static Vector3[] hcFilter(Vector3[] sv, Vector3[] pv, int[] t, float alpha, float beta)
	{
		Vector3[] wv = new Vector3[sv.Length];
		Vector3[] bv = new Vector3[sv.Length];



		// Perform Laplacian Smooth
		wv = laplacianFilter(sv, t);

		// Compute Differences
		for (int i = 0; i < wv.Length; i++)
		{
			bv[i].x = wv[i].x - (alpha * sv[i].x + (1 - alpha) * sv[i].x);
			bv[i].y = wv[i].y - (alpha * sv[i].y + (1 - alpha) * sv[i].y);
			bv[i].z = wv[i].z - (alpha * sv[i].z + (1 - alpha) * sv[i].z);
		}

		List<int> adjacentIndexes = new List<int>();

		float dx = 0.0f;
		float dy = 0.0f;
		float dz = 0.0f;

		for (int j = 0; j < bv.Length; j++)
		{
			adjacentIndexes.Clear();

			// Find the bv neighboring vertices
			adjacentIndexes = MeshUtils.findAdjacentNeighborIndexes(sv, t, sv[j]);

			dx = 0.0f;
			dy = 0.0f;
			dz = 0.0f;

			for (int k = 0; k < adjacentIndexes.Count; k++)
			{
				dx += bv[adjacentIndexes[k]].x;
				dy += bv[adjacentIndexes[k]].y;
				dz += bv[adjacentIndexes[k]].z;

			}

			wv[j].x -= beta * bv[j].x + ((1 - beta) / adjacentIndexes.Count) * dx;
			wv[j].y -= beta * bv[j].y + ((1 - beta) / adjacentIndexes.Count) * dy;
			wv[j].z -= beta * bv[j].z + ((1 - beta) / adjacentIndexes.Count) * dz;
		}

		return wv;
	}
}

public class MeshSmoothener
{
    public enum Filter { Laplacian = 1, HC = 2 };

	public static Mesh SmoothMesh(Mesh mesh, int power, Filter filterType)
    {
		for(int i = 0; i < power; ++i)
		{
			if (filterType == Filter.HC)
				mesh.vertices = SmoothFilter.hcFilter(mesh.vertices, mesh.vertices, mesh.triangles, 0.0f, 0.5f);
			if (filterType == Filter.Laplacian)
				mesh.vertices = SmoothFilter.laplacianFilter(mesh.vertices, mesh.triangles);
		}
		return mesh;
    }
}