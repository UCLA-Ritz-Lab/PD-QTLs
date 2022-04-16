#include<iostream>
#include<sstream>
using namespace std;

int datarows;
int datacols;

int main(int argc,char * argv[]){
  if (argc<3){
    cerr<<"Usage: transpose_float <datarows> <datacols>\n";
    exit(1);
  }
  datarows = atoi(argv[1]);
  datacols = atoi(argv[2]);
  cerr<<"Allocating "<<datarows<<" rows and "<<datacols<<" cols.\n";
  string *  colnames = new string[datacols];
  string *  rownames = new string[datarows];
  string line;
  getline(cin,line);
  string topleft;
  istringstream iss(line);
  iss>>topleft;
  for(int i=0;i<datacols;++i){
    iss>>colnames[i];
  }
  int rownum = 0;
  float * datamat = new float[datarows * datacols];
  while (cin){
    getline(cin,line);
    istringstream iss(line);
    iss>>rownames[rownum];
    for(int colnum = 0;colnum<datacols;++colnum){
      iss>>datamat[rownum*datacols+colnum];
    }
    cerr<<"Line "<<rownum<<" read.\n";
    ++rownum;
  }
  cerr<<"Done reading\n";
  cout<<topleft;
  for(int rownum=0;rownum<datarows;++rownum){
    cout<<"\t"<<rownames[rownum];
  }
  cout<<endl;
  for(int colnum=0;colnum<datacols;++colnum){
    cout<<colnames[colnum];
    for(int rownum=0;rownum<datarows;++rownum){
      cout<<"\t"<<datamat[rownum*datacols+colnum];
    }
    cout<<endl;
    cerr<<"Line "<<colnum<<" written.\n";
  }
  delete datamat;
  cerr<<"Done!\n";
  exit(0);
}
