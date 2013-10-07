#include "mesh.h"
#include "matrix.h"
#include <cstring>
#include <iostream>
#include <strstream>
#include <fstream>
#include <cmath>
#include <float.h>
using namespace std;

Weight currentWeight = Uniform; // Used to distinguish uniformed Laplacian & contangent weighted Laplacian.

ofstream output; // Used to output internal results for testing

/////////////////////////////////////////
// helping inline functions

inline double Cot(const Vector3d & p1, const Vector3d & p2, const Vector3d & p3) {
	Vector3d v1 = p1 - p2;
	Vector3d v2 = p3 - p2;

	v1 /= v1.L2Norm();
	v2 /= v2.L2Norm();
	double tmp = v1.Dot(v2);
	return 1.0 / tan(acos(tmp));
}

inline double Area(const Vector3d & p1, const Vector3d & p2, const Vector3d & p3) {
	Vector3d v1 = p2 - p1;
	Vector3d v2 = p3 - p1;
	return v1.Cross(v2).L2Norm() / 2.0;
}


/////////////////////////////////////////
// implementation of OneRingHEdge class
OneRingHEdge::OneRingHEdge(const Vertex * v) {
	if (v == NULL) start = next = NULL;
	else start = next = v->HalfEdge();
}

HEdge * OneRingHEdge::NextHEdge() {
	HEdge *ret = next;
	if (next && next->Prev()->Twin() != start)
		next = next->Prev()->Twin();
	else
		next = NULL;
	return ret;
}

/////////////////////////////////////////
// implementation of Mesh class
//
// function AddFace
// it's only for loading obj model, you do not need to understand it
void Mesh::AddFace(int v1, int v2, int v3) {
	int i;
	HEdge *he[3], *bhe[3];
	Vertex *v[3];
	Face *f;

	// obtain objects
	for (i=0; i<3; i++) he[i] = new HEdge();
	for (i=0; i<3; i++) bhe[i] = new HEdge(true);
	v[0] = vList[v1];
	v[1] = vList[v2];
	v[2] = vList[v3];
	f = new Face();

	// connect prev-next pointers
	SetPrevNext(he[0], he[1]);
	SetPrevNext(he[1], he[2]);
	SetPrevNext(he[2], he[0]);
	SetPrevNext(bhe[0], bhe[1]);
	SetPrevNext(bhe[1], bhe[2]);
	SetPrevNext(bhe[2], bhe[0]);

	// connect twin pointers
	SetTwin(he[0], bhe[0]);
	SetTwin(he[1], bhe[2]);
	SetTwin(he[2], bhe[1]);

	// connect start pointers for bhe
	bhe[0]->SetStart(v[1]);
	bhe[1]->SetStart(v[0]);
	bhe[2]->SetStart(v[2]);
	for (i=0; i<3; i++) he[i]->SetStart(v[i]);

	// connect start pointers
	// connect face-hedge pointers
	for (i=0; i<3; i++) {
		v[i]->SetHalfEdge(he[i]);
		v[i]->adjHEdges.push_back(he[i]);
		SetFace(f, he[i]);
	}
	v[0]->adjHEdges.push_back(bhe[1]);
	v[1]->adjHEdges.push_back(bhe[0]);
	v[2]->adjHEdges.push_back(bhe[2]);

	// mearge boundary if in need
	for (i=0; i<3; i++) {
		Vertex *start = bhe[i]->Start();
		Vertex *end   = bhe[i]->End();
		for (size_t j=0; j<end->adjHEdges.size(); j++) {
			HEdge *curr = end->adjHEdges[j];
			if (curr->IsBoundary() && curr->End()==start) {
				SetPrevNext(bhe[i]->Prev(), curr->Next());
				SetPrevNext(curr->Prev(), bhe[i]->Next());
				SetTwin(bhe[i]->Twin(), curr->Twin());
				bhe[i]->SetStart(NULL);	// mark as unused
				curr->SetStart(NULL);	// mark as unused
				break;
			}
		}
	}

	// finally add hedges and faces to list
	for (i=0; i<3; i++) heList.push_back(he[i]);
	for (i=0; i<3; i++) bheList.push_back(bhe[i]);
	fList.push_back(f);
}

// function LoadObjFile
// it's only for loading obj model, you do not need to understand it
bool Mesh::LoadObjFile(const char *filename) {
	if (filename==NULL || strlen(filename)==0) return false;
	ifstream ifs(filename);
	if (ifs.fail()) return false;

	Clear();

	char buf[1024], type[1024];
	do {
		ifs.getline(buf, 1024);
		istrstream iss(buf);
		iss >> type;

		// vertex
		if (strcmp(type, "v") == 0) {
			double x, y, z;
			iss >> x >> y >> z;
            AddVertex(new Vertex(x,y,z));
		}
		// face
		else if (strcmp(type, "f") == 0) {
			int index[3];
			iss >> index[0] >> index[1] >> index[2];
			AddFace(index[0]-1, index[1]-1, index[2]-1);
		}
	} while (!ifs.eof());
	ifs.close();

	size_t i;
	Vector3d box = this->MaxCoord() - this->MinCoord();
	for (i=0; i<vList.size(); i++) vList[i]->SetPosition(vList[i]->Position() / box.X());

	Vector3d tot;
	for (i=0; i<vList.size(); i++) tot += vList[i]->Position();
	Vector3d avg = tot / vList.size();
	for (i=0; i<vList.size(); i++) vList[i]->SetPosition(vList[i]->Position() - avg);

	HEdgeList list;
	for (i=0; i<bheList.size(); i++)
		if (bheList[i]->Start()) list.push_back(bheList[i]);
	bheList = list;

	for (i=0; i<vList.size(); i++) 
	{
		vList[i]->adjHEdges.clear();
		vList[i]->SetIndex((int)i);
		vList[i]->SetFlag(0);
	}

	return true;
}

void Mesh::DisplayMeshInfo()
{
	cout<<"/****************************/\n"<<
		  "/* Display Mesh Information */\n"<<
		  "/****************************/\n";

	cout<<"vertices # = "<<this->vList.size()<<endl;
	cout<<"non boundary half-edges # = "<<this->heList.size()<<endl;
	cout<<"boundary half-edges # = "<<this->bheList.size()<<endl;
	cout<<"faces # = "<<this->fList.size()<<endl;

	// clear bheList flag
	for(size_t i=0; i<this->bheList.size(); i++){this->bheList[i]->SetFlag(false);}

	// clear heList flag
	for(size_t i=0; i<this->heList.size(); i++){this->heList[i]->SetFlag(false);}
	
	// boundary loops
	int bloop_count = 0;
	for(size_t i=0; i<this->bheList.size(); i++){
		HEdge* curE = this->bheList[i];
		if(!curE->Flag()){// not visited loop
			bloop_count++;
			while(!curE->Flag()){ // not visited edge
				curE->SetFlag(true);
				curE->Twin()->SetFlag(true);
				curE = curE->Next();
				
				/*
				OneRingHEdge ring(curE->End());
				HEdge *curr = NULL;
				while (curr=ring.NextHEdge()) { 
					if (curr->IsBoundary() && !(curr->Flag())){
						curE=curr;
						break;
					}
				}*/
			}
		}
	}
	cout<<"boundary loop # =" << bloop_count<<endl;
	
	// clear bheList flag
	for(size_t i=0; i<this->bheList.size(); i++){this->bheList[i]->SetFlag(false);}

	// clear heList flag
	for(size_t i=0; i<this->heList.size(); i++){this->heList[i]->SetFlag(false);}

	// connected components
	int conn_count = 0;
	HEdgeList stack;
	for(size_t i=0; i<this->heList.size(); i++){
		HEdge* startE = this->heList[i];
		if(!startE->Flag()){// not visited connected component
			stack.push_back(startE);
			conn_count++;

			while(!stack.empty()){// not visited edge
				HEdge* curE = stack.back();
				curE->SetFlag(true);
				curE->Twin()->SetFlag(true);
				stack.pop_back();
			
				OneRingHEdge ring(curE->End());
				HEdge *curr = NULL;
				while (curr=ring.NextHEdge()) { 
					if ( !(curr->Flag())){
						stack.push_back(curr);
					}
				}
			}
		}
	}
	// clear bheList flag
	//for(size_t i=0; i<this->bheList.size(); i++){this->bheList[i]->SetFlag(false);}

	// clear heList flag
	//for(size_t i=0; i<this->heList.size(); i++){this->heList[i]->SetFlag(false);}

	//for(size_t i=0; i<this->vList.size(); i++){this->vList[i]->SetFlag(false);}
	//// connected components
	//int conn_count = 0;
	//

	////double test_v_count = 0;
	//for(size_t i=0; i<this->vList.size(); i++){this->vList[i]->SetFlag(false);}
	//VertexList stack;
	//for(size_t i=0; i<this->vList.size(); i++){
	//	Vertex* startV = this->vList[i];
	//	if(!startV->Flag()){
	//		//startV->SetFlag(true);
	//		stack.push_back (startV);
	//		conn_count++;

	//		while(!stack.empty()){
	//			Vertex* curV =stack.back();
	//			curV->SetFlag(true);
	//			//test_v_count+=1;
	//			//cout<<test_v_count/vList.size()<<endl;
	//			stack.pop_back();

	//			OneRingVertex ring(curV);
	//			Vertex* rV = NULL;
	//			while(rV = ring.NextVertex()){
	//				if(! rV->Flag()){
	//					//rV->SetFlag(true);
	//					stack.push_back(rV);
	//				}
	//			}
	//		}
	//	}
	//}


	cout<<"connected components # =" << conn_count<<endl;

	//for(size_t i=0; i<this->vList.size(); i++){this->vList[i]->SetFlag(false);}
	// clear bheList flag
	for(size_t i=0; i<this->bheList.size(); i++){this->bheList[i]->SetFlag(false);}
	// clear heList flag
	for(size_t i=0; i<this->heList.size(); i++){this->heList[i]->SetFlag(false);}

	// genus : v-e+f-h = 2 (c-g)
	int v = this->vList.size();
	int e = (heList.size()+bheList.size())/2.0;
	int f = this->fList.size();
	int h = bloop_count;
	int c = conn_count;

	double genus = c- (v-e+f+h)/2.0;
	cout<<"genus # = "<<genus<<endl;

}

void Mesh::HSVtoRGB( double *r, double *g, double *b, double h, double s , double v )
{
	int i;
	double f, p, q, t;
	if( s == 0 ) {
		// achromatic (grey)
		*r = *g = *b = v;
		return;
	}
	h /= 60;			// sector 0 to 5
	i = floor( h );
	f = h - i;			// factorial part of h
	p = v * ( 1 - s );
	q = v * ( 1 - s * f );
	t = v * ( 1 - s * ( 1 - f ) );
	switch( i ) {
		case 0:
			*r = v;
			*g = t;
			*b = p;
			break;
		case 1:
			*r = q;
			*g = v;
			*b = p;
			break;
		case 2:
			*r = p;
			*g = v;
			*b = t;
			break;
		case 3:
			*r = p;
			*g = q;
			*b = v;
			break;
		case 4:
			*r = t;
			*g = p;
			*b = v;
			break;
		default:		// case 5:
			*r = v;
			*g = p;
			*b = q;
			break;
	}
}
// -------------------------------------------------------
// Implement the following functions
// -------------------------------------------------------


void Mesh::ComputeVertexNormals() 
{
	// Since vertices can be classified into internal or boundary, 
	// we need to compute each class respectively.

	// 1 Compute the normals of internal vertices
	for(size_t i=0; i<this->heList.size(); i++) {
		
		Vertex* currentVertex  = heList[i]->Start();
		// Use the geometry information to compute t_1 and t_2;
		Vector3d t_1(0.0, 0.0, 0.0);
		Vector3d t_2(0.0, 0.0, 0.0);
		int k = currentVertex->Valence();
		if( k < 2) {
			cout << "Vertex Valence Computed Error!" << endl;
			return;
		}
		OneRingVertex ring(currentVertex);
		Vertex* temp = NULL;
		for(size_t i=0; i<k; i++) {
			temp = ring.NextVertex();
			t_1 = t_1 + cos(2*PI*i/k)*temp->Position();
			t_2 = t_2 + sin(2*PI*i/k)*temp->Position();
		}
		Vector3d normal = t_1.Cross(t_2);
		normal /= normal.L2Norm();
		currentVertex->SetNormal(normal);	
	}

	// 2 Compute the normals of boundary vertices
	for(size_t i=0; i<this->bheList.size(); i++) {

		Vertex* currentVertex = bheList[i]->Start();
		Vector3d t_along(0.0, 0.0, 0.0);
		Vector3d t_across(0.0, 0.0, 0.0);
		OneRingVertex ring(currentVertex);
		Vertex *p0, *p1, *pend;

		int k = currentVertex->Valence();
		if(k = 2) {
			p0 = ring.NextVertex();
			p1 = ring.NextVertex();
			t_along = p0->Position() - p1->Position();
			t_across = p0->Position() + p1->Position() - 2*currentVertex->Position();
		} else if(k = 3) {
			p0 = ring.NextVertex();
			p1 = ring.NextVertex();
			pend = ring.NextVertex();
			t_along = p0->Position() - pend->Position();
			t_across = p1->Position() - currentVertex->Position();
		} else if(k >= 4) {
			double theta = PI/(k-1);
			p0 = ring.NextVertex();
			Vertex* temp;
			for(size_t i=1; i<k-1; i++) {
				temp = ring.NextVertex();
				t_across = t_across + sin(i*theta)*temp->Position();
			}
			pend = ring.NextVertex();
			t_across = t_across*2*(cos(theta) - 1);
			t_across = t_across + sin(theta)*(p0->Position() + pend->Position());
			t_along = p0->Position() - pend->Position();
		} else {
			cout << "Vertex Valence Computed Error!" << endl;
			return;
		}
		Vector3d normal = t_along.Cross(t_across);
		normal /= normal.L2Norm();
		currentVertex->SetNormal(normal);
	}
}

void Mesh::UmbrellaSmooth() 
{
//output.open("Lc.txt");
	// We use either uniform Laplacian L_u or cotangent Laplacian L_c to do the explict smoothing.
	double lambda = 0.5;
	if(currentWeight == Uniform) {
		// We do the explict smoothing using L_u
		Vector3d *Lu = new Vector3d[vList.size()];
		for(size_t i=0; i<vList.size(); i++) {
			Vertex* v = vList[i];
			int k = v->Valence();
			OneRingVertex iterator(v);
			Lu[i] = Vector3d(0.0, 0.0, 0.0);
			Vertex *temp;
			for(size_t j=0; j<k; j++) {
				temp = iterator.NextVertex();
				Lu[i] += temp->Position() - v->Position();
			}
			Lu[i] = Lu[i]/k;
		}
		// After computing L_u for each vertex, we do the smoothing
		for(size_t i=0; i<vList.size(); i++) {
			Vector3d temp(0.0, 0.0, 0.0);
			Vertex *curr = vList[i];
			temp = curr->Position() + lambda*Lu[i];
			curr->SetPosition(temp);
		}
		cout << "Umbrella smoothing using uniform Laplacian finished.." << endl;
	} else {
		// The following code is quite similar to computing mean curvatures
		Vector3d *Lc = new Vector3d[vList.size()];
		// We do the explict smoothing using L_c
		for(size_t i=0; i<vList.size(); i++) {
			Vertex* v = vList[i];
			int k = v->Valence();
			Lc[i] = Vector3d(0.0, 0.0, 0.0);
			double cos_alpha = 0.0, cos_beta = 0.0, cot_alpha = 0.0, cot_beta = 0.0;
			double sum_w = 0.0;	// sum_w is the sum of all the weights of points shared vertex p.
			// We compute the area of triangle by using Heron's formula
			Vertex* p_pre, *p_j, *p_nex, *p_0, *p_1;
			OneRingVertex ring(v);
			p_0 = ring.NextVertex();
			p_j = p_0;
			p_1 = ring.NextVertex();
			p_nex = p_1;
			for(int j=1; j<k-1; j++) {
				p_pre = p_j;
				p_j = p_nex;
				p_nex = ring.NextVertex();

				cos_alpha = (v->Position()-p_nex->Position()).Dot(p_j->Position()-p_nex->Position()) / ((v->Position()-p_nex->Position()).L2Norm()*(p_j->Position()-p_nex->Position()).L2Norm());
				cos_beta = (v->Position()-p_pre->Position()).Dot(p_j->Position()-p_pre->Position()) / ((v->Position()-p_pre->Position()).L2Norm()*(p_j->Position()-p_pre->Position()).L2Norm());
				cot_alpha = cos_alpha / sqrt(1-cos_alpha*cos_alpha);
				cot_beta = cos_beta / sqrt(1 - cos_beta*cos_beta);
				sum_w += cot_alpha + cot_beta;
				Lc[i] += (cot_alpha+cot_beta)*(p_j->Position() - v->Position());
			}
			// Add up the calculation of start and end triangles
			p_pre = p_j;
			p_j = p_nex;
			p_nex = p_0;
			cos_alpha = (v->Position()-p_nex->Position()).Dot(p_j->Position()-p_nex->Position()) / ((v->Position()-p_nex->Position()).L2Norm()*(p_j->Position()-p_nex->Position()).L2Norm());
			cos_beta = (v->Position()-p_pre->Position()).Dot(p_j->Position()-p_pre->Position()) / ((v->Position()-p_pre->Position()).L2Norm()*(p_j->Position()-p_pre->Position()).L2Norm());
			cot_alpha = cos_alpha / sqrt(1-cos_alpha*cos_alpha);
			cot_beta = cos_beta / sqrt(1 - cos_beta*cos_beta);
			sum_w += cot_alpha+cot_beta;
			Lc[i] += (cot_alpha+cot_beta)*(p_j->Position() - v->Position());
			p_pre = p_j;
			p_j = p_nex;
			p_nex = p_1;
			cos_alpha = (v->Position()-p_nex->Position()).Dot(p_j->Position()-p_nex->Position()) / ((v->Position()-p_nex->Position()).L2Norm()*(p_j->Position()-p_nex->Position()).L2Norm());
			cos_beta = (v->Position()-p_pre->Position()).Dot(p_j->Position()-p_pre->Position()) / ((v->Position()-p_pre->Position()).L2Norm()*(p_j->Position()-p_pre->Position()).L2Norm());
			cot_alpha = cos_alpha / sqrt(1-cos_alpha*cos_alpha);
			cot_beta = cos_beta / sqrt(1 - cos_beta*cos_beta);
			sum_w += cot_alpha+cot_beta;
			Lc[i] += (cot_alpha+cot_beta)*(p_j->Position() - v->Position());
			Lc[i] = Lc[i]/sum_w;
		}
		// After computing L_c for each vertex, we do the smoothing
		for(size_t i=0; i<vList.size(); i++) {
			Vector3d temp(0.0, 0.0, 0.0);
			Vertex *curr = vList[i];
			temp = curr->Position() + lambda*Lc[i];
			curr->SetPosition(temp);
		}
		cout << "Umbrella smoothing using cotangent weight Laplacian finished.." << endl;
	}
//output.close();
}

void Mesh::ImplicitUmbrellaSmooth()
{
	double lambda = 1.0;
	Matrix *A = new Matrix(vList.size(), vList.size());
	if(currentWeight == Uniform) {
		// We use uniform Laplacian operator in this case
		// First we construct matrix A = I-\lambda*L
		for(size_t i=0; i<vList.size(); i++) {
			for(size_t j=0; j<vList.size(); j++) {
				if(i==j) {
					A->AddElement(i, j, 1-lambda);
				} else if(vList[i]->isIncidentTo(vList[j])) {
					A->AddElement(i, j, -1/vList[i]->Valence());
				} else {
					A->AddElement(i, j, 0.0);
				}
			}
		}
		// Then we perform Conjugated Gradient Method for each dimension
		
		
	} else {
		// We use cotangent weighted Laplacian operator in this case

	}
}
void Mesh::ComputeVertexCurvatures()
{
output.open("curvatures.txt");
	// 1 Compute the mean curvature of each vertex
	curvatures = new double[vList.size()];
	for(size_t i=0; i<vList.size(); i++) {

		Vertex* p = vList[i];
		int k = p->Valence();
		if(k<=2) {
			cout << "Mean Curvature Computed Error!" << endl;
			return;
		}
		Vector3d mean_curvature(0.0, 0.0, 0.0);
		double cos_alpha = 0.0, cos_beta = 0.0, cot_alpha = 0.0, cot_beta = 0.0;
		double A = 0.0;	// A is the sum of all the triangle areas shared vertex p.
		// We compute the area of triangle by using Heron's formula
		double s = 0.0;	// S = \frac{1}{2}(a + b + c)
		double a= 0.0, b = 0.0, c = 0.0; // Three edge lengths of the triangle.
		Vertex* p_pre, *p_j, *p_nex, *p_0, *p_1;
		OneRingVertex ring(p);
		p_0 = ring.NextVertex();
		p_j = p_0;
		p_1 = ring.NextVertex();
		p_nex = p_1;

		for(int j=1; j<k-1; j++) {
			p_pre = p_j;
			p_j = p_nex;
			p_nex = ring.NextVertex();

			// Calculate the interior triangle areas
			a = (p->Position() - p_j->Position()).L2Norm();
			b = (p_j->Position() - p_nex->Position()).L2Norm();
			c = (p_nex->Position() - p->Position()).L2Norm();
			s = (a + b + c)/2;
			A += sqrt(s*(s-a)*(s-b)*(s-c));

			// Calculate the curvature part
			cos_alpha = (p->Position()-p_nex->Position()).Dot(p_j->Position()-p_nex->Position()) / ((p->Position()-p_nex->Position()).L2Norm()*(p_j->Position()-p_nex->Position()).L2Norm());
			cos_beta = (p->Position()-p_pre->Position()).Dot(p_j->Position()-p_pre->Position()) / ((p->Position()-p_pre->Position()).L2Norm()*(p_j->Position()-p_pre->Position()).L2Norm());
			cot_alpha = cos_alpha / sqrt(1-cos_alpha*cos_alpha);
			cot_beta = cos_beta / sqrt(1 - cos_beta*cos_beta);
			mean_curvature += (cot_alpha+cot_beta)*(p_j->Position() - p->Position());
		}
		// Add up the calculation of start and end triangles
		p_pre = p_j;
		p_j = p_nex;
		p_nex = p_0;
		a = (p->Position() - p_j->Position()).L2Norm();
		b = (p_j->Position() - p_nex->Position()).L2Norm();
		c = (p_nex->Position() - p->Position()).L2Norm();
		s = (a + b + c)/2;
		A += sqrt(s*(s-a)*(s-b)*(s-c));
		cos_alpha = (p->Position()-p_nex->Position()).Dot(p_j->Position()-p_nex->Position()) / ((p->Position()-p_nex->Position()).L2Norm()*(p_j->Position()-p_nex->Position()).L2Norm());
		cos_beta = (p->Position()-p_pre->Position()).Dot(p_j->Position()-p_pre->Position()) / ((p->Position()-p_pre->Position()).L2Norm()*(p_j->Position()-p_pre->Position()).L2Norm());
		cot_alpha = cos_alpha / sqrt(1-cos_alpha*cos_alpha);
		cot_beta = cos_beta / sqrt(1 - cos_beta*cos_beta);
		mean_curvature += (cot_alpha+cot_beta)*(p_j->Position() - p->Position());
		p_pre = p_j;
		p_j = p_nex;
		p_nex = p_1;
		a = (p->Position() - p_j->Position()).L2Norm();
		b = (p_j->Position() - p_nex->Position()).L2Norm();
		c = (p_nex->Position() - p->Position()).L2Norm();
		s = (a + b + c)/2;
		A += sqrt(s*(s-a)*(s-b)*(s-c));
		cos_alpha = (p->Position()-p_nex->Position()).Dot(p_j->Position()-p_nex->Position()) / ((p->Position()-p_nex->Position()).L2Norm()*(p_j->Position()-p_nex->Position()).L2Norm());
		cos_beta = (p->Position()-p_pre->Position()).Dot(p_j->Position()-p_pre->Position()) / ((p->Position()-p_pre->Position()).L2Norm()*(p_j->Position()-p_pre->Position()).L2Norm());
		cot_alpha = cos_alpha / sqrt(1-cos_alpha*cos_alpha);
		cot_beta = cos_beta / sqrt(1 - cos_beta*cos_beta);
		mean_curvature += (cot_alpha+cot_beta)*(p_j->Position() - p->Position());

		// Do the division
		mean_curvature = mean_curvature/(-4*A);
		// Then we get the final mean curvature of vertex p.
		curvatures[i] = mean_curvature.L2Norm();
//output << "Curvature[" << i << "]  " << mean_curvature << endl;
	}

	// 2 Visualize the curvature by using HSV color space.
	// First we map the value of curvature to 0%-100%
	double max_curva =  0.0;
	for(size_t i=0; i<vList.size(); i++) {
		if(max_curva < curvatures[i]) max_curva = curvatures[i];
	}
	for(size_t i=0; i<vList.size(); i++) {
		curvatures[i] = curvatures[i]/max_curva;
		double red = 0.0, green = 0.0, blue = 0.0;
		double s = 1.0, v = 1.0;
		HSVtoRGB(&red, &green, &blue, curvatures[i]*360, s, v);
		//cout << "s" << curvatures[i] << endl;
		//cout << "r: " << red << ", g: " << green << ", b: " << blue << endl; 
		vList[i]->SetColor(Vector3d(red, green, blue));
	}
output.close();
}

