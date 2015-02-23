using System;
using System.Collections;
using System.Collections.Generic;

class GradientColour
{
	
	class ColourPoint           // Internal class used to store colours at different points in the gradient
	{
        public float r, g, b;   // Red, green and blue values of our colour.
	    public float val;       // Position of our colour along the gradient (between 0 and 1)

        public ColourPoint(float red, float green, float blue, float value)
        {
            r = red;
            g = green;
            b = blue;
            val = value;
        }

	}
	
	private List<ColourPoint> colour;      // An array of colour points in ascending value
	

	void Start()
	{
		createDefaultHeatMapGradient();
	}
	
	// Inserts a new colour point into its correct position:
	public void addColourPoint(float red, float green, float blue, float value)
	{
        int i = 0;
		foreach(ColourPoint col in colour)
		{
			if( value < col.val)
			{
                colour.Insert(i, new ColourPoint(red, green, blue, value));
				return;
			}
            i++;
		}
		//colour.push_back(ColourPoint(red, green,blue, value));
	}
	
	public void clearGradient()
	{
		colour.Clear();
	}
	
	// Places a 5 colour heatmap gradient into the "colour" vector:
	public void createDefaultHeatMapGradient()
	{
		colour.Clear();
		colour.Add(new ColourPoint(0, 0, 1, 0.0f));      // Blue
        colour.Add(new ColourPoint(0, 1, 1, 0.25f));     // Cyan
        colour.Add(new ColourPoint(0, 1, 0, 0.5f));      // Green
        colour.Add(new ColourPoint(1, 1, 0, 0.75f));     // Yellow
        colour.Add(new ColourPoint(1, 0, 0, 1.0f));      // Red
	}
	
	// Inputs a (value) between 0 and 1 and outputs the (red), (green) and (blue)
	//  values representing that position in the gradient.
	public void getcolourAtValue(float value, ref float red, ref float green, ref float blue)
	{
		if( colour.Count == 0)
			return;
		
        ColourPoint currC;
        ColourPoint prevC;

        for (int i = 0; i < colour.Count; i++)
		{
            currC = colour[i];
			
			if( value < currC.val)
			{
				prevC = colour[ Math.Max(0, i - 1) ];

				float valueDiff    = (prevC.val - currC.val);
				float fractBetween = (valueDiff == 0) ? 0 : (value - currC.val) / valueDiff;

				red   = (prevC.r - currC.r) * fractBetween + currC.r;
				green = (prevC.g - currC.g) * fractBetween + currC.g;
				blue  = (prevC.b - currC.b) * fractBetween + currC.b;

				return;
			}
		}

        ColourPoint last = colour[colour.Count - 1];

        red = last.r;
        green = last.g;
        blue = last.b;
		
		return;
	}

}