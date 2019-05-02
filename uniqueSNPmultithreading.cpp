#include <string>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <vector>
#include <pthread.h>
#include <mutex>
#include <math.h>
#include <thread>

using namespace std;
mutex mylock;

string getHeader(string headerLine) //just get rid of >
{
	string f;
	istringstream headerStream(headerLine);
	getline(headerStream,f,'>');
	getline(headerStream,f,'>');

return f;
}

/*
void appendSequence(vec, string header, string sequence)
{
	auto fai=fa.find(header);
	if(fai==fa.end())
		fa[header]=sequence;
	else
		(fai->second)+=sequence;
}*/
	
struct UniqueSNP
{
public:
	string strain;
	int position;
	char type;
};

void uSNPThread(vector<string> header, vector<string> sequence, int starti, int endi, ofstream &out, string runMode,vector<int> bucket) //bucket need to have 5 elements
{
	for(int i=starti; i<endi; ++i)
	{
		//five buckets, occupied by first occurence; if conflicted, cover bucket
		bucket[0]=-1;
		bucket[1]=-1;
		bucket[2]=-1;
		bucket[3]=-1;
		bucket[4]=-1;	

		for(int seqi=0; seqi<header.size(); seqi++)
		{
			
			if(sequence[seqi][i]=='A' || sequence[seqi][i]=='a')
			{
				if(bucket[0]==-1)
					bucket[0]=seqi;
				else
					bucket[0]=-2;
			}	
			else if(sequence[seqi][i]=='T' || sequence[seqi][i]=='t')
			{
				if(bucket[1]==-1)
					bucket[1]=seqi;
				else
					bucket[1]=-2;
			}	
			else if(sequence[seqi][i]=='G' || sequence[seqi][i]=='g')
			{
				if(bucket[2]==-1)
					bucket[2]=seqi;
				else
					bucket[2]=-2;
			}	
			else if(sequence[seqi][i]=='C' || sequence[seqi][i]=='c')
			{
				if(bucket[3]==-1)
					bucket[3]=seqi;
				else
					bucket[3]=-2;
			}	
			else if(sequence[seqi][i]=='-')
			{
				if(bucket[4]==-1)
					bucket[4]=seqi;
				else
					bucket[4]=-2;
			}	
		}
		
		if(runMode=="2" && bucket[4]!=-1) //if runmode2, skip gaps
			continue;
		mylock.lock();	
		if(bucket[0]>=0)
		{
			out<<header[bucket[0]]<<"\t"<<i+1 <<"\tA"<<endl;
		}
		if(bucket[1]>=0)
		{
			out<<header[bucket[1]]<<"\t"<<i+1 <<"\tT"<<endl;
		}
		if(bucket[2]>=0)
		{
			out<<header[bucket[2]]<<"\t"<<i+1 <<"\tG"<<endl;
		}
		if(bucket[3]>=0)
		{
			out<<header[bucket[3]]<<"\t"<<i+1 <<"\tC"<<endl;
		}
		if(runMode=="0")
		{
			if(bucket[4]>=0)
			{
				out<<header[bucket[4]]<<"\t"<<i+1 <<"\t-"<<endl;
			}
		}
		mylock.unlock();

	}
}	 

int main(int argc, char *argv[])
{
	string alignmentName=argv[1];
	string runMode=argv[2];
	string threadm=argv[3];
	int p=stoi(threadm);


	//read fasta store into map
	cout<<"\t1-Reading in fasta..."<<endl;
	ifstream alignment(alignmentName.c_str());
	string line;
	vector<string> header;
	vector<string> sequence;
	//map<string, string> fa; //strain - sequence
	int seqi=-1;

	for(;getline(alignment,line);)
	{
		if(line=="") continue;
		if(line[0]=='>')
		{
			header.push_back(getHeader(line));
			sequence.push_back("");
			seqi++;
		}
		else
		{
			sequence[seqi]+=line;
		}
	}

	//record unique snps
	cout<<"\t2-Recording unique snvs..."<<endl;
	//vector<UniqueSNP> uniqueSNP;
	ofstream out("uniqueSNPs");
	out<<"strain\tposition\tbasetype"<<endl;

	vector<int> bucket={-1,-1,-1,-1,-1};
	int numThread=p;
	if(sequence[0].size() < p)
		numThread=sequence[0].size();

	
	int len=ceil((double)sequence[0].size()/(double)numThread);

	vector<vector<int>> buckets;
	for(int bs=0; bs<numThread;bs++)
	{
		buckets.push_back(bucket);
	}

	thread t[numThread];
	for(int threadi=0; threadi<numThread; threadi++)
	{
		int starti=threadi*len;
		int endi=(threadi+1)*len-1;
		if(endi>=sequence[0].size())
			endi=sequence[0].size()-1;

		t[threadi]=thread(uSNPThread,header,sequence, starti, endi, std::ref(out), runMode,buckets[threadi]);
		//linei++;
	}

	for(int threadi=0; threadi<numThread; threadi++)
	{
		t[threadi].join();
	}	
	

/*
	for(int i=0; i<((fa.begin())->second).size(); ++i)
	{
		//five buckets, occupied by first occurence; if conflicted, cover bucket
		string bucketA="";
		string bucketT="";
		string bucketG="";
		string bucketC="";
		string bucketGap="";

		for(auto seq : fa)
		{
			if(seq.second[i]=='A' || seq.second[i]=='a')
			{
				if(bucketA=="")
					bucketA=seq.first;
				else
					bucketA="coveredBucket";
			}	
			else if(seq.second[i]=='T' || seq.second[i]=='t')
			{
				if(bucketT=="")
					bucketT=seq.first;
				else
					bucketT="coveredBucket";
			}
			else if(seq.second[i]=='G' || seq.second[i]=='g')
			{
				if(bucketG=="")
					bucketG=seq.first;
				else
					bucketG="coveredBucket";
			}
			else if(seq.second[i]=='C' || seq.second[i]=='c')
			{
				if(bucketC=="")
					bucketC=seq.first;
				else
					bucketC="coveredBucket";
			}
			else if(seq.second[i]=='-')
			{				
				if(bucketGap=="")
					bucketGap=seq.first;
				else
					bucketGap="coveredBucket";
			}
		}

		if(runMode=="2" && bucketGap!="") //if runmode2, skip gaps
			continue;
		
		if(bucketA!="" && bucketA!="coveredBucket")
		{
			UniqueSNP newUniqueSNP;
			newUniqueSNP.strain=bucketA;
			newUniqueSNP.position=i+1;
			newUniqueSNP.type='A';
			uniqueSNP.push_back(newUniqueSNP);
		}
		if(bucketT!="" && bucketT!="coveredBucket")
		{
			UniqueSNP newUniqueSNP;
			newUniqueSNP.strain=bucketT;
			newUniqueSNP.position=i+1;
			newUniqueSNP.type='T';
			uniqueSNP.push_back(newUniqueSNP);
		}
		if(bucketG!="" && bucketG!="coveredBucket")
		{
			UniqueSNP newUniqueSNP;
			newUniqueSNP.strain=bucketG;
			newUniqueSNP.position=i+1;
			newUniqueSNP.type='G';
			uniqueSNP.push_back(newUniqueSNP);
		}
		if(bucketC!="" && bucketC!="coveredBucket")
		{
			UniqueSNP newUniqueSNP;
			newUniqueSNP.strain=bucketC;
			newUniqueSNP.position=i+1;
			newUniqueSNP.type='C';
			uniqueSNP.push_back(newUniqueSNP);
		}

		if(runMode=="0")
		{
			if(bucketGap!="" && bucketGap!="coveredBucket")
			{
				UniqueSNP newUniqueSNP;
				newUniqueSNP.strain=bucketGap;
				newUniqueSNP.position=i+1;
				newUniqueSNP.type='-';
				uniqueSNP.push_back(newUniqueSNP);
			}
		}

	}
*/
/*
	for(auto us : uniqueSNP)
	{
		output<<us.strain<<"\t"<<us.position<<"\t"<<us.type<<endl;
	}
*/
}
