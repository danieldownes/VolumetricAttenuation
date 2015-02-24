// Converted from UnityScript to C# at http://www.M2H.nl/files/js_to_c.php - by Mike Hergaarden
// Do test the code! You usually need to change a few small bits.

using UnityEngine;
using System.Collections;

public class AttenuationRays : MonoBehaviour
{
	public GameObject lightmeshholder;

	public float gradientBand;
	//gradientBandUpper / Lower
	
	private int RaysToShoot=512; //64; 128; 1024; 
	private int distance=15;
	private Vector3[] vertices;
	private Vector2[] vertices2d;
	private int[] triangles;
	//private Vector3[] vertices2;
	private Mesh mesh;


	public GradientColour gradient;
	

	void Start ()
	{
		
		//vertices = new Vector3[RaysToShoot];
		vertices2d = new Vector2[RaysToShoot];
		//triangles = new int[RaysToShoot];
		//	vertices2 = new Vector3[4];
		mesh= lightmeshholder.GetComponent<MeshFilter>().mesh;

		gradient = new GradientColour();
		gradient.init();
		gradient.createDefaultHeatMapGradient();
	}
	
	void Update()
	{
		vertices = mesh.vertices;
		
		float angle = 0;

		RaycastHit hit;
		RaycastHit hit2;
		RaycastHit hitB;

		float rayStrength;

		// Shoot rays
		for (int i=0; i < RaysToShoot;i++)
		{
			float x = Mathf.Sin(angle);
			float y = Mathf.Cos(angle);

			rayStrength = 1;
			
			angle += 2 * Mathf.PI / RaysToShoot;
			
			Vector3 dir = new Vector3(x,0,y);

			// First Hit
			if (Physics.Raycast(transform.position, dir, out hit, distance))
			{
				Debug.DrawLine (transform.position, hit.point, new Color(1,0,0,0.7f));


				// Check for object behind object of interest

				Vector3 hp = hit.point;
				hp += dir * 0.01f;
				if (Physics.Raycast(hp, dir, out hit2, distance)) 
				{
					Vector3 hp2 = hit2.point;
					//Debug.DrawLine(hp, hp2, new Color(1,0.5f,0,0.7f));
					
					// Trace back to object of interest
					//dir = -dir; // * 0.01f;
					if (Physics.Raycast(hp2, -dir, out hitB, distance)) 
					{
						Vector3 hpb = hitB.point;
						//hitB.point
						//Debug.DrawLine (hp2, hp2 + (Vector3.forward * 0.3f), Color(1,1,1,0.7f));
						//Debug.DrawLine (hpb, hpb + (Vector3.left * 0.3f), Color(1,1,1,0.7f));

						float dist  = Vector3.Distance(hit.point, hpb);

						rayStrength -= dist * gradientBand;

						//hp2 = Mathf.Max(0f, 1.0f);

						Color col = new Color(0, 0, 0, 1);

						gradient.getcolourAtValue(rayStrength, ref col.r, ref col.g, ref col.b);

						Debug.DrawLine (hp2, hpb, col);
					}
					
				}
				
				//FIXME_VAR_TYPE tmp= lightmeshholder.transform.InverseTransformPoint(hit.point);
				//vertices[i] = Vector3(tmp.x,0.5f,tmp.z);
			}
			else
			{ 
				// no hit
				//Debug.DrawRay (transform.position, dir*distance, Color(1,1,0,1));
				//FIXME_VAR_TYPE tmp2= lightmeshholder.transform.InverseTransformPoint(lightmeshholder.transform.position+dir);
				//vertices[i] = Vector3(tmp2.x,0.5f,tmp2.z);
			}
		}
		
		// last vertice is at the player location (center point)
		//vertices[i] = lightmeshholder.transform.InverseTransformPoint(transform.position);
		
		//mesh.vertices = vertices;
	}
	
	
	
	
	void  BuildMesh (){
		//TO CONSIDER
		// dont cast if not moved?
		// build prelook-array of hit points/pixels/areas?
		// skip duplicate hit points (compare previous)
		// always same amount of vertices, no need create new mesh?..but need to triangulate or not??
		
		float angle = 0;
		for (int i=0;i < RaysToShoot;i++)
		{
			float  x= Mathf.Sin(angle);
			float  y= Mathf.Cos(angle);
			angle += 2 * Mathf.PI / RaysToShoot;
			
			Vector3 dir = new Vector3(x, 0, y);
			RaycastHit hit;
			if (Physics.Raycast (transform.position, dir, out hit, distance)) 
			{
				Debug.DrawLine (transform.position, hit.point, new Color(1,1,0,1));
				
				//// bounce ray, ignore last hit object?
				//// reflect more than 1 ray?
				//Vector3 dir2 = Vector3.Reflect(dir, hit.normal);
				//RaycastHit hit2;
				//if (Physics.Raycast (hit.point, dir2, hit2, distance/4)) 
				//{
				//Debug.DrawLine (hit.point, hit2.point, Color(1,0,0,1));
				//}else{ // we might not hit anything, because of short bounce distance
				//	Debug.DrawRay (hit.point, dir2*(distance/4), Color(1,0,0,1));
				//}
				
				
				//Vector2 tmp= lightmeshholder.transform.InverseTransformPoint(hit.point);
				//vertices2d[i] = Vector2(tmp.x,tmp.z);
				
			}
			else
			{ 
				// no hit
				//			Debug.DrawRay (transform.position, dir*distance, Color(1,1,0,1));
				
				//FIXME_VAR_TYPE tmp2= lightmeshholder.transform.InverseTransformPoint(lightmeshholder.transform.position+dir);
				//vertices2d[i] = Vector2(tmp2.x,tmp2.z);
				
			}
		}
		
		
		/*
	// triangulate.cs
//    Triangulator tr = new Triangulator(vertices2d);
//    int[] indices = tr.Triangulate();
	
	// build mesh
    Vector2[] uvs = new Vector2[vertices2d.Length+1];
    Vector3[] newvertices = new Vector3[vertices2d.Length+1];
    for (int n = 0; n<newvertices.Length-1;n++) 
	{
        newvertices[n] = new Vector3(vertices2d[n].x, 0, vertices2d[n].y);

	// create some uv's for the mesh?
	// uvs[n] = vertices2d[n];
		
    }
    
	//print("len"+newvertices.Length+" n:"+n);
	
	triangles = new int[newvertices.Length*3];
	    
	// triangle list
	i = -1;
	for (n=0;n<triangles.length-3;n+=3)
	{
		i++;
		triangles[n] = newvertices.Length-1;
		if (i>=newvertices.Length)
		{
			triangles[n+1] = 0;
			//print ("hit:"+i);
		}else{
			triangles[n+1] = i+1;
		}
		triangles[n+2] = i;
	}    
    i++;
    
	// central point
	newvertices[newvertices.Length-1] = new Vector3(0,0,0);
	triangles[triangles.length-3] = newvertices.Length-1;
	triangles[triangles.length-2] = 0;
	triangles[triangles.length-1] = i-1;
   
    // Create the mesh
    //Mesh msh = new Mesh();
    mesh.vertices = newvertices;
    //DD/mesh.triangles = triangles;
    //DD/mesh.uv = uvs;
	
//    mesh.RecalculateNormals(); // need?
//    mesh.RecalculateBounds(); // need ?

	// last triangles
//	triangles[i+1] = 0;
//	triangles[i+2] = 0;
//	triangles[i+1] = 0;

	//triangles.Reverse();

//	mesh.vertices = newvertices;
//	mesh.triangles = triangles;

	// not every frame? clear texture before take new shot?
//	if (grab>10) GrabToTexture();
//	grab++;
	*/
	}
	
}