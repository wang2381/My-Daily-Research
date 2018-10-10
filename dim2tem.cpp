#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
 
using namespace std;
 
typedef vector<float>          OneDimension;
typedef vector<OneDimension> TwoDimension;
 


void Display2 (TwoDimension const& t)
{
//   float sum[22359]={0};
/*
  
   for (TwoDimension::size_type i = 0; i < t.size(); ++i) //row 488160
   {
      for (OneDimension::size_type j = 0; j < t[i].size(); ++j) //column 720
      {
        sum[j]=sum[j]+t[i][j];
      }
   }
*/
   FILE *fd;
   fd=fopen("TS14K148.564.5after.txt","w");
   float temp, sum=0, average, max, min;

   for (int j=0;j<1997;j++)
   {
       fprintf(fd,"%.1f,%.1f, %s ,%d,%d,",t[0][0],t[0][1],"Ts",1294,1900+j);
       max=t[j][3];
       min=t[j][3];
       for (int k=0;k<12;k++)
       {
           sum+=t[j][k+3];
           if (t[j][k+3]>max){max=t[j][k+3];}
           if (t[j][k+3]<min){min=t[j][k+3];}
           
       }
       average=sum/12;
       fprintf(fd,"%.1f,%.1f,%.2f,%.1f,",sum,max,average,min);

       sum=0;            
       for (int k=0;k<12;k++)
       {
           if(k<11)
           {fprintf(fd,"%.1f,",t[j][k+3]);}
           else{fprintf(fd,"%.1f, %s\n",t[j][k+3],"Globe");
       }
   }
}
}

int main (int const argc, char const* argv[])
{
   TwoDimension v2;
   ifstream     f("TS_14k_148.5_64.5.txt");
 
   if (!f)
   {
      cerr << "open a.txt failed" << endl;
      return (-1);
   }
 
   string s;
   float    t;
   while (getline(f, s))
   {
      istringstream is(s);
      OneDimension  v1;
 
      while (is >> t)
      {
         v1.push_back(t);
      }
      v2.push_back(v1);
   }
   f.close();
   cout<<v2[1][1]<<endl;
   Display2(v2);
 
   return (0);
}
