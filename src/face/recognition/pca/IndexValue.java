package face.recognition.pca;
/*-----------------------------------------------------------------------------
File         : IndexValue
 Author      : Engin DURMAZ
 Date        : 23.04.2019
 Description : Class to store a label and its index together
-----------------------------------------------------------------------------*/
public class IndexValue implements Comparable
{
	int index;
	double value;

	public IndexValue(int i, double v) 
	{
		index = i;
		value = v;
	}

	public int compareTo(Object obj) 
	{
		if(obj instanceof IndexValue)
		{
			double target = ((IndexValue) obj).value;
			
			if (value > target)
			{
				return -1;
			}
			else if (value < target)
			{
				return 1;
			}
			return 0;
		}

		return -1000;
	}
}
