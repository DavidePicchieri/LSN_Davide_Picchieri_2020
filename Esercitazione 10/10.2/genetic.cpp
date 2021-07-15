#include "genetic.h"

using namespace std;



Specimen :: Specimen(int ncities, Random * rnd){
  for(int i=0; i<ncities; i++) 
    m_cities.push_back(i);
  m_rnd=rnd;
}

Specimen :: ~Specimen(){}


Specimen& Specimen :: operator=( Specimen& a) // Overload
{
  for(int i=0; i<m_cities.size(); i++){
    m_cities[i]=a.GetCities()[i];
  }
  m_distance = a.m_distance;
  m_rnd = a.m_rnd;
  m_gen = a.m_gen;
  
  return *this;
}

void Specimen :: Check (){
  vector<int> sorted;
  for(int i=0; i<m_cities.size(); i++) sorted.push_back(m_cities[i]);
  sort(sorted.begin(), sorted.end());
  for(int i=0; i< m_cities.size()-1; i++){
    if(sorted[i]==sorted[i+1]){
      cout << "Two cities are the same!" << endl;
      abort();
    }
  }
}


void Specimen :: Print(string filename){
  ofstream Output;
  Output.open(filename, ios::app);
  for(int i=0; i<m_cities.size(); i++){
    Output << m_cities[i] << endl;
  }
  Output << endl << endl;
}

int Specimen :: Pbc(int i){ //Periodic boundary conditions without the 0th element
  if(i >  m_cities.size()){
    i = i - m_cities.size()+1; 
  }
  else if(i == m_cities.size()){
    i = i - m_cities.size()+1;
  } 
  else if(i == 0){
    i = i + m_cities.size()-1;
  } 
  else if(i < 0){  
    i = i + m_cities.size()-1;
  }
  return i;
}

void Specimen :: Swap(int i){
  if(i>=m_cities.size()) {
    cout << "Error i >= N" << endl;
    abort();
  }
  
  if(i==0) {
    cout << "Error i=0" << endl;
    abort();
  }
  int appo=m_cities[i];
  m_cities[i]=m_cities[Pbc(i+1)];
  m_cities[Pbc(i+1)]=appo;
}

void Specimen :: Shift(int i, int groupsize, int shift){
  if(groupsize>=m_cities.size()-1) {
    cout << "Error groupsize>=N-1" << endl;
    abort();
  }
  if(i>=m_cities.size() or i==0) {
    cout << "Error i>=N or i=0" << endl;
    abort();
  }
  if(shift>=m_cities.size()-1) {
    cout << "Error shift>=N-1" << endl;
    abort();
  }
  
  if(shift+groupsize>=m_cities.size()-1) {
    cout << "Error shift+groupsize>=N" << endl;
    abort();
  }
  
  int appo1[shift+groupsize];
  int appo2[shift+groupsize];
  
  for(int j=0; j<groupsize; j++) appo1[j]=m_cities[Pbc(i+j)];
  for(int j=0; j<shift; j++) appo2[j]=m_cities[Pbc(Pbc(i+groupsize)+j)];
  for(int j=0; j<shift; j++) m_cities[Pbc(i+j)]=appo2[j];
  for(int j=0; j<groupsize; j++) m_cities[Pbc(Pbc(i+shift)+j)]=appo1[j];
  
}

void Specimen :: SwapGroup(int i, int groupsize){
  if(groupsize>=m_cities.size()/2) {
    cout << "Error groupsize>=N/2" << endl;
    abort();
  }
  if(i>=m_cities.size() or i==0) {
    cout << "Error i>=N or i=0" << endl;
    abort();
  }
  
  int appo1[groupsize];
  int appo2[groupsize];
  
  for(int j=0; j<groupsize; j++){
    appo1[j]=m_cities[Pbc(i+j)];
    appo2[j]=m_cities[Pbc(Pbc(i+groupsize)+j)];
    m_cities[Pbc(Pbc(i+groupsize)+j)]=appo1[j];
    m_cities[Pbc(i+j)]=appo2[j];
  }
}

void Specimen :: Inversion(int i, int groupsize){
  if(groupsize>=m_cities.size()-1) {
    cout << "Error groupsize>=N-1" << endl;
    abort();
  }
  if(i>=m_cities.size() or i==0) {
    cout << "Error i>=N or i=0" << endl;
    abort();
  }
  
  vector<int> appo;
  for(int j=0; j<groupsize; j++) appo.push_back(m_cities[Pbc(i+j)]);
  std::reverse(appo.begin(), appo.end());
  for(int j=0; j<groupsize; j++) m_cities[Pbc(i+j)]=appo[j];
}

///////////////////SPECIES////////////////////////////////

Species :: Species(int popsize, Random * rnd, int ncities){
    m_popsize = popsize;
    m_rnd = rnd;
    m_ncities = ncities;
    for(int i=0; i<m_popsize; i++){
        Specimen appo(m_ncities, m_rnd);
        m_pop.push_back(appo);
        m_popnew.push_back(appo);
    }
    m_ngen = 0;
}
Species :: ~Species(){}

/////////////////////SPECIES : PRINT/////////////////////////

void Species :: Print(int index, string filename){
  m_pop[index].Print(filename);
}
void Species :: PrintBest(string filename){
  ofstream Output;
  Output.open(filename);  
  for(int i=0; i<m_ncities; i++){
    Output << m_besttrip[i] << "\t" << m_pos[0][m_besttrip[i]] << "\t" << m_pos[1][m_besttrip[i]] << endl;
  }
  Output << m_besttrip.size() << "\t" << m_pos[0][0] << "\t" << m_pos[1][0] << endl;
  Output.close();
}

void Species :: PrintDistance(string filename){ //Overloaded to print whole population or single individual
  ofstream Output;
  Output.open(filename);
  if (!Output.is_open()){
      abort;
  }
  for (int i=0; i<m_popsize; i++){
      Output << m_pop[i].GetDistance() << endl;
  }
  Output.close();
}
void Species :: PrintDistance(int index, string filename){
  ofstream Output;
  Output.open(filename, ios::app);
  Output << m_ngen << setw(12) << m_pop[index].GetDistance() << endl;
  Output.close();
}
void Species :: PrintAve(int nspecimen, string filename){
  double ave = 0;
  if (nspecimen > m_popsize){
    cout << "Insert < " << m_popsize << " as argument of the function PrintAve" << endl;
    abort;
  }
  for (int i=0; i<nspecimen; i++){
    ave += m_pop[i].GetDistance();
  }
  ave /= double(nspecimen);
  ofstream Output;
  Output.open(filename, ios::app);
  Output << m_ngen << "\t" << ave << endl;
  Output.close();
}

void Species :: PrintBestSA(string filename){
  ofstream Output;
  Output.open(filename);
  
  for(int i=0; i<m_ncities; i++){
    Output << m_pop[0].GetCities()[i] << setw(12) << m_pos[0][m_pop[0].GetCities()[i]] << setw(12) << m_pos[1][m_pop[0].GetCities()[i]] << endl;
  }
  Output << m_pop[0].GetCities().size() << setw(12) << m_pos[0][0] << setw(12) << m_pos[1][0] << endl;
}


void Species :: PrintBestRank(string filename){
  ofstream Output;
  Output.open(filename);  
  for(int i=0; i<m_ncities; i++){
    Output << m_besttrip[i] << "\t" << m_pos[0][m_besttrip[i]] << "\t" << m_pos[1][m_besttrip[i]] << endl;
  }
  Output << m_besttrip.size() << "\t" << m_pos[0][0] << "\t" << m_pos[1][0] << endl;F
}


/////////////////////////////SPECIES : METHODS///////////////////////////
/////////////////////////////SPECIES : METHODS///////////////////////////
/////////////////////////////SPECIES : METHODS///////////////////////////


void Species::Copy(int i, int j){
    m_pop[i]=m_pop[j];
}  


///////////////////SPECIES : INITIALIZE POSITION////////////////////////
void Species :: SetPosCircle(){
  vector<double> x, y;
  double angle, r=1;
  x.push_back(0); //initial position always set in (0,1)
  y.push_back(1);

  for (int i=1; i<m_ncities; i++){
      angle = m_rnd->Rannyu(0,2*M_PI);
      x.push_back(r*cos(angle));
      y.push_back(r*sin(angle));
  }
    
  m_pos.push_back(x);
  m_pos.push_back(y);


  for(int i=0; i<1000; i++){ //Swapping many times the cities in the different specimen to create the initial population
    for(int j=0; j<m_popsize; j++){
      m_pop[j].Swap(int(m_rnd->Rannyu(1, m_ncities)));
    }
  }
  m_besttrip = m_pop[0].GetCities();

}

void Species :: SetPosSquare(){
  double side = 1;
  vector<double> x, y;
  for (int i=0; i<m_ncities; i++){
      x.push_back(m_rnd->Rannyu(-side, side));
      y.push_back(m_rnd->Rannyu(-side, side));
  }

  m_pos.push_back(x);
  m_pos.push_back(y);

  for(int i=0; i<1000; i++){ //Swapping many times the cities in the different specimen to create the initial population
    for(int j=0; j<m_popsize; j++){
      m_pop[j].Swap(int(m_rnd->Rannyu(1, m_ncities)));
    }
  }
  m_besttrip = m_pop[0].GetCities();
}

//////////////////////////////////SPECIES : SET & GET////////////////////////////////

void Species :: SetCoeff(double coeff){
    m_ann_coeff = coeff;
}
void Species :: SetTemp(double temp){
  m_temp = temp;
}
void Species :: SetPos(int i, int j, double pos){
  m_pos[i][j]=pos;
}

double Species :: GetAlpha(){
  return fmin(1., exp(-1./m_temp*(m_pop[1].GetDistance()-m_pop[0].GetDistance())));
}

//////////////////////////////////SPECIES : COMPLEX METHODS////////////////////////////////

void Species :: Measure(){
  int appo1, appo2;
  double distance, dx, dy;

  for(int i=0; i<m_popsize; i++){
    distance=0;

    for(int j=0; j<m_ncities-1; j++){
      appo1 = m_pop[i].GetCities()[j];
      appo2 = m_pop[i].GetCities()[j+1];
      dx = m_pos[0][appo1] - m_pos[0][appo2];
      dy = m_pos[1][appo1] - m_pos[1][appo2];
      distance += sqrt(dx*dx+dy*dy);
    }

    appo1 = m_pop[i].GetCities()[m_ncities-1];
    appo2 = m_pop[i].GetCities()[0];
    dx = m_pos[0][appo1] - m_pos[0][appo2];
    dy = m_pos[1][appo1] - m_pos[1][appo2];
    distance += sqrt(dx*dx+dy*dy);
    
    m_pop[i].SetDistance(distance);
  }

}
void Species :: Sort(){
  int index;
  m_besttrip = m_pop[0].GetCities(); 
  for (int i = 0; i < m_popsize-1; i++){
    index = i;
    for (int j = i+1; j < m_popsize; j++){
      if (m_pop[j].GetDistance() < m_pop[index].GetDistance()) 
        index = j;
    }
    Specimen appo = m_pop[i];
    m_pop[i] = m_pop[index];
    m_pop[index] = appo;
  }
  if(m_pop[0].GetDistance()<m_mindist){
    m_besttrip = m_pop[0].GetCities();
    m_mindist = m_pop[0].GetDistance();
    m_bestgen = m_ngen;
  }
}
void Species :: Selection(double powerlaw){
  int index1, index2;
  for(int i=0; i<m_popsize; i=i+2){
    
    index1=int(m_popsize*pow(m_rnd->Rannyu(),powerlaw));
    do {
      index2=int(m_popsize*pow(m_rnd->Rannyu(),powerlaw));
    }
    while(index2==index1);
    
    m_popnew[i]=m_pop[index1];
    m_popnew[i+1]=m_pop[index2];
  }
}
void Species :: Crossover(int cut, double chance){
  int count1, count2;
  for(int i=0; i<m_popsize; i=i+2){
      m_pop[i]=m_popnew[i];
      m_pop[i+1]=m_popnew[i+1];
      if(m_rnd->Rannyu() < chance){
        
        int missing1[m_ncities-cut];
        int missing2[m_ncities-cut];
        for(int j=0; j<m_ncities-cut; j++){
          missing1[j]=m_popnew[i].GetCities()[j+cut];
          missing2[j]=m_popnew[i+1].GetCities()[j+cut];
        }    
        count1=cut;
        count2=cut;
        for(int j=1; j<m_ncities; j++){ 
          for(int k=0; k<m_ncities-cut; k++){
            if(m_popnew[i+1].GetCities()[j] == missing1[k]){
              m_pop[i].SetCity(count1, missing1[k]);
              count1++;
            }
            if(m_popnew[i].GetCities()[j] == missing2[k]){
              m_pop[i+1].SetCity(count2, missing2[k]);
              count2++;
            }
          }
        }
        for(int k=0; k<m_popsize; k++) m_pop[k].Check();
      }
    }
  m_ngen=m_ngen+1;

}
void Species :: Mutation(double chance){
  int groupsize, shift;

  for(int i=0; i<m_popsize; i++){  //Swap
    if(m_rnd->Rannyu()<chance){
      m_pop[i].Swap(int(m_rnd->Rannyu(1, m_ncities)));
    }
  }

  for(int i=0; i<m_popsize; i++){  //Shift
      if(m_rnd->Rannyu()<chance){
        do{
          groupsize = int(m_rnd->Rannyu(1,m_ncities-1));
          shift = int(m_rnd->Rannyu(1,m_ncities-1));
        }while(groupsize + shift >=(m_ncities-1));

        m_pop[i].Shift(int(m_rnd->Rannyu(1, m_ncities)), groupsize, shift);
    }
   }

  for(int i=0; i<m_popsize; i++){  //SwapGroup
      if(m_rnd->Rannyu()<chance){
      groupsize = int(m_rnd->Rannyu(1, m_ncities/2));
      m_pop[i].SwapGroup(int(m_rnd->Rannyu(1, m_ncities)), groupsize);
    }
  }
  
  for(int i=0; i<m_popsize; i++){  //Inversion
      if(m_rnd->Rannyu()<chance){
      groupsize = int(m_rnd->Rannyu(1, m_ncities-1));
      m_pop[i].Inversion(int(m_rnd->Rannyu(1, m_ncities)), groupsize);
    }
  }

  for(int i=0; i<m_popsize; i++) m_pop[i].Check();
}


void Species :: Annealing(){
  m_temp *= m_ann_coeff;
}


void Species :: Accept(){
  if(m_rnd->Rannyu()<this->GetAlpha()){
    m_pop[0]=m_pop[1];
  }else{
    m_pop[1]=m_pop[0];
  }
  m_ngen=m_ngen+1;
}

void Species :: Rearrange(){
  int j = m_rnd->Rannyu(1, m_ncities);
  m_pop[1].Swap(j);

  j=int(m_rnd->Rannyu(1, m_ncities));
  int m = m_rnd->Rannyu(1, m_ncities-1);
  m_pop[1].Inversion(j, m);
  
  m_pop[1].Check();
}
 
/////////////////////////METHODS FOE 10.2///////////////////////////// 

void Species :: MigrationDirection(int * vec, int & target){//randomly choose which rank will talk to which one
  vector<int> appo;
  target = m_rnd->Rannyu(1,4);

  for(int i=0; i<4; i++) 
    appo.push_back(i);

  appo.erase(appo.begin()+target);
  appo.erase(appo.begin());

  vec[0]=appo[0];
  vec[1]=appo[1];
}

void Species :: Migration(int & rank, int target, int * vec){//Perform the migration
  int appo1[m_ncities], appo2[m_ncities];
  int generation, savegen;
  int itag = 1, itag2 = 2, itag3 = 3, itag4 = 4;
  
  for(int i=0; i<m_ncities; i++){ //put the best path in a temporary vector
    appo1[i]=m_pop[0].GetCities()[i];
  }

  
  MPI_Status status_1, status_2, status_3, status_4;
 
  if(rank==0){          //send the best to another rank and receive back their best
    cout << rank << " Migrating to :"  << target << endl;
    MPI_Send(appo1, 32, MPI_INTEGER, target, itag, MPI_COMM_WORLD);
    MPI_Recv(appo2, 32, MPI_INTEGER, target, itag2, MPI_COMM_WORLD, &status_2);
  }else if(rank == target){
    cout << rank << " Migrating to :"  << 0 << endl;
    MPI_Recv(appo2, 32, MPI_INTEGER, 0, itag, MPI_COMM_WORLD, &status_1);
    MPI_Send(appo1, 32, MPI_INTEGER, 0, itag2, MPI_COMM_WORLD);
  }
  
  if(rank==vec[0]){
    cout << rank << " Migrating to :"  << vec[1] << endl;
    MPI_Send(appo1, m_ncities, MPI_INTEGER, vec[1], itag3, MPI_COMM_WORLD);
    MPI_Recv (appo2, m_ncities, MPI_INTEGER, vec[1], itag4, MPI_COMM_WORLD, &status_4);
  }else if(rank==vec[1]){
    cout << rank << " Migrating to :"  << vec[0] << endl;
    MPI_Recv(appo2, m_ncities, MPI_INTEGER, vec[0], itag3, MPI_COMM_WORLD, &status_3);
    MPI_Send(appo1, m_ncities, MPI_INTEGER, vec[0], itag4, MPI_COMM_WORLD);
  }
  
  for(int i=0; i<m_ncities; i++){  //save the received path
    m_pop[0].GetCities()[i]=appo2[i];
  }
}

void Species :: Best(int & rank){ //find which irank has the bast path and share it
  
  int rece0[m_ncities];
  int gen0, gen1, gen2, gen3;
  double length1, length2, length3;
  int send1[m_ncities];
  int send2[m_ncities];
  int send3[m_ncities];
  int itag1=1;
  int itag2=2;
  int itag3=3;
  int itag4=4;
  int itag5=5;
  int itag6=6;
  int itag7=7;
  int itag8=8;
  int itag9=9;
  
  
  for(int i=0; i<m_ncities; i++){ //preparing to send and receive the best paths and their length 
    rece0[i]=m_besttrip[i];
  }
  gen0=m_bestgen;

  if (rank == 1){
    length1 = m_mindist;
    for(int i=0; i<m_ncities; i++){
      send1[i]=m_besttrip[i];
    }
    gen1 = m_bestgen;
  }
  else if (rank == 2){
    length2 = m_mindist;
    for(int i=0; i<m_ncities; i++){
      send2[i]=m_besttrip[i];
    }
    gen2 = m_bestgen;
  }
  else if (rank == 3){
    length3 = m_mindist;
    for(int i=0; i<m_ncities; i++){
      send3[i]=m_besttrip[i];
    }    
    gen3 = m_bestgen;
  }
   //sending and receiving the best paths and their length
    
  MPI_Status stat1, stat2, stat3, stat4, stat5, stat6, stat7, stat8, stat9;
  
  if(rank==1){
    MPI_Send(rece0,32,MPI_INTEGER,0,itag1,MPI_COMM_WORLD);
  }else if(rank==0){
    MPI_Recv(send1,32,MPI_INTEGER,1,itag1, MPI_COMM_WORLD,&stat1);
  }
  
  if(rank==0){
    MPI_Recv(send2,32,MPI_INTEGER,2,itag2, MPI_COMM_WORLD,&stat2);
  }else if(rank==2){
   MPI_Send(rece0,32,MPI_INTEGER,0,itag2,MPI_COMM_WORLD);
  }
  
  if(rank==0){
    MPI_Recv(send3,32,MPI_INTEGER,3,itag3, MPI_COMM_WORLD,&stat3);
  }else if(rank==3){
    MPI_Send(rece0,32,MPI_INTEGER,0,itag3,MPI_COMM_WORLD);
  }
  
  if(rank==0){
    MPI_Recv(&gen1,1,MPI_INTEGER,1,itag4, MPI_COMM_WORLD,&stat4);
    MPI_Recv(&length1,1,MPI_DOUBLE,1,itag5, MPI_COMM_WORLD,&stat5);
  }else if(rank==1){
    MPI_Send(&gen0,1,MPI_INTEGER,0,itag4,MPI_COMM_WORLD);
    MPI_Send(&m_mindist,1,MPI_DOUBLE,0,itag5,MPI_COMM_WORLD);
  }

  if(rank==0){
    MPI_Recv(&gen2,1,MPI_INTEGER,2,itag6, MPI_COMM_WORLD,&stat6);
    MPI_Recv(&length2,1,MPI_DOUBLE,2,itag7, MPI_COMM_WORLD,&stat7);
  }else if(rank==2){
    MPI_Send(&gen0,1,MPI_INTEGER,0,itag6,MPI_COMM_WORLD);
    MPI_Send(&m_mindist,1,MPI_DOUBLE,0,itag7,MPI_COMM_WORLD);

  }
  
  if(rank==0){
    MPI_Recv(&gen3,1,MPI_INTEGER,3,itag8, MPI_COMM_WORLD,&stat8);
    MPI_Recv(&length3,1,MPI_DOUBLE,3,itag9, MPI_COMM_WORLD,&stat9);
  }else if(rank==3){
    MPI_Send(&gen0,1,MPI_INTEGER,0,itag8,MPI_COMM_WORLD);
    MPI_Send(&m_mindist,1,MPI_DOUBLE,0,itag9,MPI_COMM_WORLD);
  }
  double mindist;
  int best_generation; //cheching which one is the best path among them all 
  if(rank==0){
    if(fmin(fmin(fmin(m_mindist,length1),length2),length3)==length1){
      for(int i=0; i<m_ncities; i++){
        m_besttrip[i]=send1[i];
      }
      mindist=length1;
      best_generation=gen1;
    }
    if(fmin(fmin(fmin(m_mindist,length1),length2),length3)==length2){
      for(int i=0; i<m_ncities; i++){
        m_besttrip[i]=send2[i];
      }
      mindist=length2;
      best_generation=gen2;
    }
    if(fmin(fmin(fmin(m_mindist,length1),length2),length3)==length3){
      for(int i=0; i<m_ncities; i++){
        m_besttrip[i]=send3[i];
      }
      mindist=length3;
      best_generation=gen3;
    }
    if(fmin(fmin(fmin(m_mindist,length1),length2),length3)==m_mindist){
      mindist=m_mindist;
      best_generation=m_bestgen;
    }
  }
  
}
