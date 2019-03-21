#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define OK 666
#define ERR 999
////////////////////////////
//////////数据结构//////////
////////////////////////////
//DBG图节点
typedef struct DBGNode
{
	int geneId;
	struct EdgeNode *link;
	int linkNum;
}DBG;

//节点地址
typedef struct IndexNode
{
	struct DBGNode *address;
}Index;

//记录相连节点的地址和权重
typedef struct EdgeNode
{
	int geneId;
	int weight;
	struct DBGNode *node;
	struct EdgeNode *next;
}Edge;

//Edge的地址列表
typedef struct EdgeListNode
{
	struct EdgeNode* node;
	struct EdgeListNode* next;
}EdgeList;

////////////////////////////
///////////子函数///////////
////////////////////////////
//参数交互
void error()
{
	fprintf(stderr,"Trace back to the chromosome ancestor\n");
	fprintf(stderr,"Author:chenhaixin@genomics.cn\n");
	fprintf(stderr,"Version:1.0.0\n");
	fprintf(stderr,"Usage:tb2ca -x 1000 -f inputfile -m 1\n ");
	fprintf(stderr,"\t-f\tinput file.\n");
	fprintf(stderr,"\t-x\tmax id.\n");
	fprintf(stderr,"\t-m\tpick up mothod[default 1].\n");
	fprintf(stderr,"\t\t\t1.don't pick up.\n");
	fprintf(stderr,"\t\t\t2.ramdon.\n");
	fprintf(stderr,"\t\t\t3.other.\n");
	fprintf(stderr,"Good Luck.\n");
	exit(0);
}

//判断该gene ID是否已经生成节点，返回节点地址
//如果没有，新建gene 节点
DBG *getNodeAddress(Index *index,int geneID)
{
	if(index[geneID].address)
	{
		return index[geneID].address;
	}
	else
	{
		DBG *newNode = (DBG *)malloc(sizeof(DBG));
		newNode->geneId = geneID;
		newNode->link = NULL;
		newNode->linkNum = 0;
		index[geneID].address = newNode;
		return newNode;
	}
}

//找到该gene的对应的权重节点,返回地址
//如果没有，需要新建
void addEdgeWeight(DBG *startNode,DBG *endNode,int geneId,int weight)
{
	//edge 头节点
	//temp 当前节点
	//last 上一个节点 
	Edge * last = startNode->link;
	Edge * temp = last;
	while(temp != NULL)
	{
		if(temp->geneId == geneId)
		{
			temp->weight += weight;
			return;
		}
		last = temp;
		temp = temp->next;
	}
	
	//新建节点
	Edge *newNode = (Edge *)malloc(sizeof(Edge));
	newNode->geneId = geneId;
	newNode->weight = weight;
	newNode->node = endNode;
	newNode->next = NULL;

	if(last == NULL)
	{
		//只有一个节点
		startNode->link = newNode;
	}
	else
	{
		//还有其它节点，查到链表末尾
		last->next = newNode;
	}
	startNode->linkNum += 1;//新增了一个节点
	return;
}

//添加两个节点的权重,如果没有，需要新建
void addWeight(int lastGeneId,int tempGeneId,DBG *lastAddress,DBG *tempAddress,int weight)
{
	addEdgeWeight(lastAddress,tempAddress,tempGeneId,weight);
	addEdgeWeight(tempAddress,lastAddress,lastGeneId,weight);
}
//构图
int drawDBG(Index *index,char file[],int maxGeneID)
{
	FILE *fp;
	if((fp = fopen(file,"r")) == NULL)
	{
		fprintf(stderr,"%s file open error. please check!\n",file);
		return ERR;
	}
	while(!feof(fp))
	{
		int weight;
		int lastGeneId;
		int tempGeneId;
		DBG *lastAddress = NULL;
		DBG *tempAddress = NULL;

		char c;
		fscanf(fp,"%d",&weight);
		c=fgetc(fp);
		fscanf(fp,"%d",&lastGeneId);
		c=fgetc(fp);
		if(feof(fp)){fprintf(stdout,"finish reads file.\n");break;}//feof不能很好地判断文件是否结束，还要在判断一次

		lastAddress = getNodeAddress(index,lastGeneId);
		while(c != '\n' && !feof(fp))
		{
			//生成节点和添加权重信息
			fscanf(fp,"%d",&tempGeneId);
			c=fgetc(fp);
			tempAddress = getNodeAddress(index,tempGeneId);
			addWeight(lastGeneId,tempGeneId,lastAddress,tempAddress,weight);
			lastGeneId = tempGeneId;
			lastAddress = tempAddress;
			//printf("%d%c",lastGeneId,c);
		}
		//换行
		//fprintf(stdout,"weight=%d\n",weight);
		//fprintf(stdout,"finish one line.\n");
	}
	fclose(fp);
	return OK;
}

//随机遍历一条染色体
void traverse(Index *index,DBG *head)
{
	if(head->linkNum ==0)
	{
		return;
	}

	fprintf(stdout,"geneId1=%d,linkNum=%d,weight=%d,geneId2=%d\n",head->geneId,head->linkNum,head->link->weight,head->link->geneId);
	head->linkNum = 0;//不再输出,告诉本次输出
	index[head->link->geneId].address = NULL;//不再输出,告诉主函数

	//找到尚未输出的,
	Edge *tempE = head->link;
	while(tempE && tempE->node && tempE->node->linkNum == 0)
	{
		tempE = tempE->next;
	}
	if(tempE && tempE->node && tempE->node->linkNum !=0){traverse(index,tempE->node);}
	return;
}

//删除边的动作
//需要删除该节点和指向节点的Edge节点
void deleteEdgeNode(DBG *node,int geneId)
{
	node->linkNum--;//减少一条边
	Edge *tempEdge = node->link;

	if(tempEdge->geneId == geneId)//删掉头节点
	{
		//fprintf(stdout,"remove,%d,%d,%d\n",node->geneId,tempEdge->weight,geneId);
		node->link = tempEdge->next;
		tempEdge->next = NULL;
		tempEdge->node = NULL;
		free(tempEdge);
		return;
	}
	else
	{
		Edge *lastEdge = tempEdge;
		tempEdge = tempEdge->next;
		//找到定位
		while(tempEdge->geneId != geneId)
		{
			lastEdge = tempEdge;
			tempEdge = tempEdge->next;
		}
		//fprintf(stdout,"remove,%d,%d,%d\n",node->geneId,tempEdge->weight,geneId);
		lastEdge->next = tempEdge->next;
		tempEdge->next = NULL;
		tempEdge->node = NULL;
		free(tempEdge);
		return;
	}
}

//批量删除动作
void deleteEdgeList(EdgeList *head,DBG *node)
{
	//先删除指向节点的权重边
	EdgeList *deleteNode = head->next;
	while(deleteNode)
	{
		//指向节点的node，出发节点的id
		deleteEdgeNode(deleteNode->node->node,node->geneId);
		deleteNode = deleteNode->next;
	}
	//然后删除出发节点的权重边
	deleteNode = head->next;
	while(deleteNode)
	{
		//出发节点的node，指向节点的id
		deleteEdgeNode(node,deleteNode->node->geneId);
		deleteNode = deleteNode->next;
	}
	return;
}

//在待删除的列表中挽留其中的1或2个
//第三个参数是挽留的个数
//mode=1 全部删掉
//mode=2 随机保留一个
//mode=3 
void pickUpEdge(EdgeList *weekHead,int mode,int luckNum)
{
	if(mode == 1)
	{
		return;
	}
	else if(mode == 2)
	{
		EdgeList *temp = weekHead->next;
		if(temp == NULL){return;}
		weekHead->next = temp->next;
		temp->node = NULL;
		temp->next = NULL;
		free(temp);
		if(luckNum == 2){pickUpEdge(weekHead,mode,1);}
		return;
	}
	else if(mode == 3)
	{
		return;
	}
	else
	{
		return;
	}
}
//遍历index，按照权重依次从低到高的顺序，迭代进行删除低权重的边，使得边的数量<=2
//1.判断边的数量，如果<=2则跳过
//2.1 如果比tempWeight大的边>=2，则删掉所有权重小的边
//2.2 如果比tempWeight大的边=1，则需要保留一条或不保留
//2.2.1 随机保留一条
//2.2.2 都删掉(结果会很碎）
//2.2.3 待定
//2.3 如果比tempWeight大的边=0，则需要保留两条或不保留
//2.3.1 随机保留一条
//2.3.2 都删掉
//2.3.3 待定
void removeWeekEdge(Index *index,int weight,int maxGeneId,int mode)
{
	int x;
	for(x=0;x<maxGeneId;x++)
	{
		DBG *node = index[x].address;
		if(node == NULL){continue;}
	
		if(node->linkNum <=2){continue;}//1.
		//找到比temp weight大的节点有几个
		Edge *tempEdge = node->link;
		EdgeList *weekHead = (EdgeList *)malloc(sizeof(EdgeList));//用于记录权重小于或等于当前权重的节点地址，头节点不写内容
		EdgeList *weekTemp = weekHead;
		weekHead->node = NULL;
		weekHead->next = NULL;
		int bigWeightEdgeNum = 0;
		int weekWeightEdgeNum = 0;
		while(tempEdge)
		{
			if(tempEdge->weight > weight)
			{
				bigWeightEdgeNum++;
			}
			else
			{
				weekWeightEdgeNum++;
				EdgeList *temp = (EdgeList *)malloc(sizeof(EdgeList));
				temp->node = tempEdge;
				temp->next = NULL;
				weekTemp->next = temp;
				weekTemp = temp;
			}
			tempEdge = tempEdge->next;
		}
		
		//如果大于或等于2，则全部删除,无需操作
		//如果等于1，则最多只能保留一个
		if(bigWeightEdgeNum == 1)
		{
			pickUpEdge(weekHead,mode,1);
		}
		//如果没有，则最多保留两个
		if(bigWeightEdgeNum == 0)
		{
			pickUpEdge(weekHead,mode,2);
		}

		deleteEdgeList(weekHead,node);
	}
}

////////////////////////////
///////////主函数///////////
////////////////////////////
int main(int argc,char *argv[])
{
	////////参数输入
	int maxGeneID = 0;
	char *inputfile = NULL;
	int mode = 1;

	int x = 0;
	for(x = 1;x < argc;x += 2)
	{
		if(strcmp(argv[x], "-f") == 0 )
		{
			inputfile = argv[x+1];
	        }
		else if(strcmp(argv[x], "-m") == 0)
		{
			mode = atoi(argv[x+1]);
		}
		else if(strcmp(argv[x], "-x") == 0)
		{
			maxGeneID = atoi(argv[x+1]);
		}
		else
		{
			error();
		}
	}
	fprintf(stdout,"maxGeneID=%d\n",maxGeneID);
	fprintf(stdout,"inputfile=%s\n",inputfile);
	maxGeneID += 1;//兼容gene id = 0；gene id就是其节点在地址数组的下标
	if(maxGeneID == 0 || inputfile == NULL){error();}
	//开始构图
	Index *index = (Index *)malloc(sizeof(Index) * (maxGeneID));
	for(x=0;x<maxGeneID;x++)
	{
		index[x].address = NULL;
	}
	if(drawDBG(index,inputfile,maxGeneID) == OK)
	{
		fprintf(stdout,"DBG draw finish.\n");
	}
	else
	{
		fprintf(stderr,"Some error happen when draw DBG\n.");
		exit(-1);
	}
	//迭代断开权重小的边
	int weight = 1;
	//traverse(index[10].address);
	removeWeekEdge(index,weight,maxGeneID,mode);

	for(x=1;x<maxGeneID;x++)
	{
		if(index[x].address)
		{
			fprintf(stdout,"###############################\n");
			traverse(index,index[x].address);
		}
	}

	return 0;
}
