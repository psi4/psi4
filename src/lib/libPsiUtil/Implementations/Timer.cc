/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */
#include<iomanip>
#include "../Timer.h"
#include "../Table.h"
#include "../../libparallel2/LibParallel2.h"
#include "psi4-dec.h"
using std::cout;
using std::endl;

namespace psi {
typedef std::pair<double,int> MyTime_t;
typedef PsiMap<TimeTypes,MyTime_t> Res_t;
typedef boost::timer::cpu_times TimeStruct_t;
typedef boost::timer::cpu_timer Timer_t;

Time_t TimeByType(const TimeTypes Type, const TimeStruct_t& TimeStruct) {
   Time_t returnvalue;
   switch (Type) {
      case (WALL): {returnvalue=TimeStruct.wall;break;}
      case (CPU): {returnvalue=TimeStruct.user;break;}
      case (SYSTEM): {returnvalue=TimeStruct.system;break;}
      default: {throw PSIEXCEPTION("Unrecognized Time Type");break;}
   }
   return returnvalue;
}

template<int NPoints>
MyTime_t DetermineRes(TimeTypes Type){
   if (!boost::chrono::high_resolution_clock::is_steady)
      (*outfile)<<"Warning!!! Your clock is not a steady clock.  "
            "Your timings may not be reliable."<<std::endl;
   TimeStruct_t start_time;
   Time_t current_time,AvgTime=0;
   start_time.clear();
   Timer_t cpu;
   for (int i=0; i<NPoints; ++i) {
      cpu.start();
      start_time=cpu.elapsed();
      current_time=TimeByType(Type,start_time);
      while(current_time==TimeByType(Type,start_time))
         current_time=TimeByType(Type,cpu.elapsed());
      AvgTime+=current_time-TimeByType(Type,start_time);
   }
   double Time=AvgTime/NPoints;
   int value=0;
   if(Time/1.0e9>1e-3)value=3;
   else if(Time/1.0e9>1e-6)value=6;
   else value=9;
   return MyTime_t(Time/1.0e9,value);
}

Res_t InitRes(){
   Res_t ReturnValue;
   ReturnValue[WALL]=DetermineRes<1>(WALL);
   ReturnValue[CPU]=DetermineRes<1>(CPU);
   ReturnValue[SYSTEM]=DetermineRes<1>(SYSTEM);
   return ReturnValue;
}

namespace TimeConsts{static  Res_t Resolutions=InitRes();}

TimeValue::TimeValue(TimeTypes Type):
      Resolution_(TimeConsts::Resolutions[Type]){}

std::string TimeValue::PrintOut(const int i)const{
   std::stringstream Result;
   Result<<std::setprecision(Time_.second)<<Time_.first<<std::fixed;
   if(i>0)Result<<" +/- "<<Resolution_.first;
   return Result.str();
}

bool SmartTimer::IsStopped() const {return Timer_.is_stopped();}

void SmartTimer::Start(){Stop();Timer_.start();StartTimes_=Timer_.elapsed();}

void SmartTimer::Resume(){Timer_.resume();}

void SmartTimer::Stop() {Timer_.stop();StopTimes_=Timer_.elapsed();}

TimeValue SmartTimer::GetTime(TimeTypes Type) const {
   TimeStruct_t temp=(IsStopped()? StopTimes_ : Timer_.elapsed());
   Time_t diff;
   diff=TimeByType(Type, temp);
   diff-=TimeByType(Type, StartTimes_);
   TimeValue Return(Type);
   int Unit=Return.Resolution_.second;
   Return.Time_.first=diff/1e9;
   Return.Time_.second=Unit;
   return Return;
}

SmartTimer::SmartTimer(const std::string& Name, bool IsParallel) :
      IsParallel_(IsParallel),Comm_(WorldComm->GetComm()){
   Name_<<Name;
   this->Start();
}

void FillArray(std::vector<double>& MyTimes,
              const SmartTimer& Timer,const TimeTypes Type){
   int offset=((int)Type)*3;
   MyTimes[offset+0]=Timer.GetTime(Type).Time_.first;
   MyTimes[offset+1]=Timer.GetTime(Type).Time_.second;
   MyTimes[offset+2]=Timer.GetTime(Type).Resolution_.first;
}
std::string SmartTimer::PrintOut() const {
   std::stringstream Result;
   const int NProc=Comm_->NProc(), DataPieces=9;
   std::vector<double> AllTimes(DataPieces*NProc), MyTimes(DataPieces);
   for(int i=0;i<3;i++)
      FillArray(MyTimes,*this,(TimeTypes)i);
   if(IsParallel_)Comm_->AllGather(&MyTimes[0],DataPieces,&AllTimes[0]);
   else AllTimes=MyTimes;
   Result<<Name_.str()<<std::endl;
   std::vector<std::string> ProcessNames;
   std::vector<std::string> FinalTimes;
   for(int i=0;i<NProc;i++){
      if(IsParallel_){
         std::stringstream tempstr;
         tempstr<<"Processor "<<i;
         ProcessNames.push_back(tempstr.str());
      }
      else ProcessNames.push_back("Total:");
      for(int j=0;j<3;j++){
         int offset=i*DataPieces+j*3;
         TimeValue temp((TimeTypes)j);
         temp.Time_.first=AllTimes[offset];
         temp.Time_.second=
               temp.Resolution_.second=(int)AllTimes[offset+1];
         temp.Resolution_.first=AllTimes[offset+2];
         FinalTimes.push_back(temp.PrintOut(j==0?0:1));
      }
   }
   if(IsParallel_){
      ProcessNames.push_back("% Imbalance :");
      for(int i=0;i<3;i++){
         double Imbalance=ComputeLoadImbalance((TimeTypes)i,AllTimes);
         std::stringstream ToString;
         ToString<<std::setprecision(2)<<std::fixed<<Imbalance;
         FinalTimes.push_back(ToString.str());
      }
   }
   typedef TableColumn<std::string> StrCol;
   StrCol RowTitles("",&ProcessNames[0],1,0,'\0','|');
   StrCol WallTimes("Wall (s)",&FinalTimes[0],3);
   StrCol CPUTimes("CPU (s)",&FinalTimes[1],3);
   StrCol SysTimes("System (s)",&FinalTimes[2],3);
   Table<StrCol,StrCol,StrCol,StrCol>
   MyTable(NProc+(IsParallel_?1:0),RowTitles,WallTimes,CPUTimes,SysTimes);
   MyTable.SetBorder(TOP,'*');
   MyTable.SetBorder(BOTTOM,'*');
   Result<<MyTable.GetTable()<<std::endl;
   return Result.str();
}

double SmartTimer::ComputeLoadImbalance(TimeTypes Type,std::vector<double>& Times)const {
   double Average=0.0,Max=-1.0;
   int NProcs=Times.size()/9;
   for(int i=0;i<NProcs;i++){
      double time=Times[i*9+(int)Type*3];
      Average+=time;
      if(Max<time)Max=time;
   }
   return (Average>0.0?(Max*NProcs/Average-1.0)*100:0.0);
}

} //End namespace
