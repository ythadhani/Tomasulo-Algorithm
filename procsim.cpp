#include "procsim.hpp"
#include <iostream>
#include <deque>
#include <algorithm>
using namespace std;

/**
 * Subroutine for initializing the processor. You many add and initialize any global or heap
 * variables as needed.
 * XXX: You're responsible for completing this routine
 *
 * @r  Number of result buses
 * @k0 Number of k0 FUs
 * @k1 Number of k1 FUs
 * @k2 Number of k2 FUs
 * @f Number of instructions to fetch
 */
struct fetched_instructions 
{
  uint64_t dest_tag;
  uint64_t address;  //decimal equivalent of address received from trace 
  uint64_t fu;
  int64_t dest_reg;
  int64_t src1_reg;
  int64_t src2_reg;
};

struct schedule_queue        //Ready bits for sources have not been used, instead tag values have been reset to zero to indicate availability
{
  uint64_t fu;
  uint64_t dest_tag;
  int64_t dest_reg;
  uint64_t src1_tag;
  int64_t src1_reg;
  uint64_t src2_tag;
  int64_t src2_reg;
  bool filled_bit=false;      // To check whether a row in the schedule queue actually contains a valid entry
  bool fired=false;	      // Fired bit is set once the instruction has both source operands available and the required FU is available
  bool completed=false;	      // Completed bit is set once the instruction has been broadcast on the result bus and can be deleted in the next cycle
};

struct functional_units
{
	bool busy_bit;
	uint64_t fu_type;
	uint64_t dest_tag;
	int64_t dest_reg;
	uint64_t time_spent=0;
};

struct result_bus
{
	bool busy_bit;
	uint64_t tag;
	int64_t dest_reg;
};

struct register_file
{
	uint64_t dest_tag;
	bool busy_bit;
};

struct output                    //In order to print cycle by cycle behavior, we store the results in array of structs which is printed once the last inst. retires
{
	uint64_t instruction_number;
	uint64_t fetched;
	uint64_t dispatched;
	uint64_t scheduled;
	uint64_t executed;
	uint64_t status_update;
};

schedule_queue* s_queue;
functional_units* f_unit;
result_bus* rb;
deque<fetched_instructions> dispatch_queue;
output* out;
uint64_t R,K0,K1,K2,F;
uint64_t tag =0;
uint64_t rows;
uint64_t total_disp_queue_size = 0;
uint64_t total_instructions_dispatched =0;
uint64_t total_instructions_fired =0;
uint64_t total_instructions_retired =0;
uint64_t maximum_dispatch_size = 0;
register_file* reg_file;
uint64_t cycle_count=0;
uint64_t number_of_retires = 0;
void fetch_to_dispatch();
void dispatch_to_schedule();
void schedule_to_fu();
void status_update();
void fu_to_resultbus();
void update_scheduling_queue();

void setup_proc(uint64_t r, uint64_t k0, uint64_t k1, uint64_t k2, uint64_t f) 
{
	R=r;
	K0=k0;
	K1=k1;
	K2=k2;
	F=f;
	rows = 2*(K0+K1+K2);
	s_queue  = new schedule_queue[rows];
	f_unit = new functional_units[K0+K1+K2];
	rb = new result_bus[R];
	reg_file = new register_file[128];
	out = new output[100008];
	for(int i=0;i<128;i++)
	{	
		reg_file[i].dest_tag = 0;	
		reg_file[i].busy_bit = 0;
	}

	for(uint64_t i=0;i<(K0+K1+K2);i++)
	{
		if(i<K0)
		f_unit[i].fu_type = 0;
		else if(i<(K0+K1))
		f_unit[i].fu_type = 1;
		else if(i<(K0+K1+K2))
		f_unit[i].fu_type = 2;
	}
}

/**
 * Subroutine that simulates the processor.
 *   The processor should fetch instructions as appropriate, until all instructions have executed
 * XXX: You're responsible for completing this routine
 *
 * @p_stats Pointer to the statistics structure
 */
void run_proc(proc_stats_t* p_stats)
{
	while(number_of_retires<100000)      //The cycle count keeps getting incremented until the last instruction has retired
	{	
		cycle_count++;
			
		fu_to_resultbus();

		schedule_to_fu();

		dispatch_to_schedule();

		update_scheduling_queue();
	
		status_update();
	
		fetch_to_dispatch();

		total_disp_queue_size = total_disp_queue_size + dispatch_queue.size();
	
		maximum_dispatch_size = max(dispatch_queue.size(),maximum_dispatch_size);
	}
}

void fetch_to_dispatch()                      //Fetches 4 instructions every cycle and populates the dispatch queue, a deque has been used to model the same
{
	for(uint64_t i=0;i<F;i++)
	{
		proc_inst_t p_inst;
		bool flag = read_instruction(&p_inst);
		if(flag==true)                // Keep reading until End of File has been reached
		{
			fetched_instructions inst;		
			tag = tag+1;
			out[tag].instruction_number = tag;
			out[tag].fetched = cycle_count;
			inst.dest_tag = tag;
			inst.address = p_inst.instruction_address;
			inst.fu = p_inst.op_code;
			if(p_inst.op_code == -1)
			{
				inst.fu = 1;
			}
			inst.dest_reg = p_inst.dest_reg;
			inst.src1_reg = p_inst.src_reg[0];
			inst.src2_reg = p_inst.src_reg[1];
			dispatch_queue.push_back(inst);
		}
		else
		{
			break;
		}
	}
		
}

void dispatch_to_schedule()                   //Transfers instructions to schedule queue if there are vacant slots 
{	
	for(uint64_t j=0;j<rows;j++)
	{
		if(s_queue[j].filled_bit==0 && !dispatch_queue.empty())
		{		
			s_queue[j].dest_tag = dispatch_queue[0].dest_tag;
			out[dispatch_queue[0].dest_tag].dispatched = cycle_count;				
			s_queue[j].fu =  dispatch_queue[0].fu;
			s_queue[j].filled_bit =1;
			s_queue[j].dest_reg = dispatch_queue[0].dest_reg;
			s_queue[j].src1_reg = dispatch_queue[0].src1_reg;
			s_queue[j].src2_reg = dispatch_queue[0].src2_reg;
			s_queue[j].completed = 0;
			s_queue[j].fired = 0;
			if(dispatch_queue[0].src1_reg!=-1 && reg_file[dispatch_queue[0].src1_reg].busy_bit ==1)
			{
				s_queue[j].src1_tag =  reg_file[dispatch_queue[0].src1_reg].dest_tag;
			}
			else
			{
				s_queue[j].src1_tag = 0;
			}
			if(dispatch_queue[0].src2_reg!=-1 && reg_file[dispatch_queue[0].src2_reg].busy_bit ==1)
			{
				s_queue[j].src2_tag =  reg_file[dispatch_queue[0].src2_reg].dest_tag;
			}
			else
			{
				s_queue[j].src2_tag = 0;
			}
			if(dispatch_queue[0].dest_reg!=-1)    //The register number is used to index to the register file array, hence -1 is not permissible
			{
				reg_file[dispatch_queue[0].dest_reg].dest_tag = dispatch_queue[0].dest_tag;
				reg_file[dispatch_queue[0].dest_reg].busy_bit = 1;
			}
			
			dispatch_queue.pop_front();
		}
	}			
}

void schedule_to_fu()   //Transfers instructions from schedule queue to FU if corresponding FU is free and source operands are available
{
	for(uint64_t i=0; i<rows; i++)
	{	
		if(s_queue[i].filled_bit ==1 && s_queue[i].fired ==0 && s_queue[i].src1_tag == 0 &&  s_queue[i].src2_tag ==0) 
		{
			for(uint64_t j=0; j<(K0+K1+K2); j++)
			{	
				if(f_unit[j].fu_type==s_queue[i].fu && f_unit[j].busy_bit==0 )
				{
					f_unit[j].busy_bit =1;
					s_queue[i].fired = 1;
					f_unit[j].dest_tag = s_queue[i].dest_tag;
					f_unit[j].dest_reg = s_queue[i].dest_reg;
					out[f_unit[j].dest_tag].scheduled = cycle_count;
					break;
				}
				
			}
		}
	}
}

void fu_to_resultbus()		//Broadcasts results of executed instructions on the result bus  
{
	for(uint64_t i=0;i<R;i++)
	{
		rb[i].busy_bit =0;
		rb[i].dest_reg =0;
		rb[i].tag =0;
		uint64_t max_time_spent;	
		uint64_t min_tag;
		int64_t min_tag_index =-1;
		for(uint64_t j=0;j<(K0+K1+K2);j++)
		{
			if(f_unit[j].busy_bit==1)
			{
				max_time_spent= f_unit[j].time_spent;
				min_tag= f_unit[j].dest_tag;
				break;
			}
		}
			
		for(uint64_t j=0;j<(K0+K1+K2);j++)  //First priority given to instructions that have spent maximum time in the FU
		{
			if(f_unit[j].busy_bit==1 && max_time_spent<=f_unit[j].time_spent)
			{
				max_time_spent = f_unit[j].time_spent;
				min_tag = f_unit[j].dest_tag;
				min_tag_index = j;
			}
		}
		for(uint64_t j=0;j<(K0+K1+K2);j++)   //Thereafter if time spent is same then follow tag order
		{
			if(f_unit[j].busy_bit==1 && f_unit[j].time_spent==max_time_spent && min_tag>=f_unit[j].dest_tag)
			{
				min_tag = f_unit[j].dest_tag;
				min_tag_index = j;
			}
		}
		if(min_tag_index !=-1)
		{
			rb[i].busy_bit =1;
			rb[i].dest_reg = f_unit[min_tag_index].dest_reg;
			rb[i].tag = f_unit[min_tag_index].dest_tag;
			f_unit[min_tag_index].dest_reg =0;
			f_unit[min_tag_index].dest_tag =0;
			f_unit[min_tag_index].busy_bit =0;
			f_unit[min_tag_index].time_spent =0;	
			if(rb[i].dest_reg!=-1 && reg_file[rb[i].dest_reg].dest_tag == rb[i].tag)
			{
				reg_file[rb[i].dest_reg].dest_tag = 0;
				reg_file[rb[i].dest_reg].busy_bit = 0;
			}
		}
	}

	for(uint64_t i=0;i<(K0+K1+K2);i++)
	{
		if(f_unit[i].busy_bit==1)
		{
			f_unit[i].time_spent = f_unit[i].time_spent+1;
		}
	}
}

void update_scheduling_queue()      //This updates the scheduling queue with the results broadcast on the result bus
{
	for(uint64_t i=0; i<rows; i++)
	{	
		for(uint64_t j=0; j<R; j++)
		{
			if(s_queue[i].filled_bit==1 && s_queue[i].src1_tag== rb[j].tag)
			{
				s_queue[i].src1_tag = 0;
			}
			if(s_queue[i].filled_bit==1 && s_queue[i].src2_tag== rb[j].tag)
			{
				s_queue[i].src2_tag = 0;
			}
		}
	}
}

void status_update()       //Deletes instructions from the scheduling queue once marked as complete
{
	uint64_t i=0;
	while(i<rows)   
	{	
		if(s_queue[i].completed == 1 && s_queue[i].filled_bit ==1)
		{	
			number_of_retires = number_of_retires + 1;
			out[s_queue[i].dest_tag].status_update = cycle_count;
			s_queue[i].fu =0;	
			s_queue[i].dest_tag = 0;
			s_queue[i].dest_reg = 0;
			s_queue[i].src1_tag = 0;
			s_queue[i].src1_reg = 0;
			s_queue[i].src2_tag = 0;
			s_queue[i].src2_reg = 0;
			s_queue[i].filled_bit = 0;
			s_queue[i].fired = 0;
			s_queue[i].completed = 0;			
			for(uint64_t j=i; j<(rows-1); j++)   //Once an inst. is deleted shift all instructions upwards to fill the slot thereby maintaining tag order
			{
				s_queue[j].fu = s_queue[j+1].fu;
				s_queue[j].dest_tag = s_queue[j+1].dest_tag;
				s_queue[j].dest_reg = s_queue[j+1].dest_reg;
				s_queue[j].src1_tag = s_queue[j+1].src1_tag;
				s_queue[j].src1_reg = s_queue[j+1].src1_reg;
				s_queue[j].src2_tag = s_queue[j+1].src2_tag;
				s_queue[j].src2_reg = s_queue[j+1].src2_reg;
				s_queue[j].filled_bit = s_queue[j+1].filled_bit;
				s_queue[j].fired= s_queue[j+1].fired;
				s_queue[j].completed= s_queue[j+1].completed;

				s_queue[j+1].fu =0;	
				s_queue[j+1].dest_tag = 0;
				s_queue[j+1].dest_reg = 0;
				s_queue[j+1].src1_tag = 0;
				s_queue[j+1].src1_reg = 0;
				s_queue[j+1].src2_tag = 0;
				s_queue[j+1].src2_reg = 0;
				s_queue[j+1].filled_bit = 0;
				s_queue[j+1].fired = 0;
				s_queue[j+1].completed = 0;
			}
		}
		else
		{
			i = i+1;
		}
	}
	
	for(uint64_t i=0; i<rows; i++)        //if instruction has been broadcast on the result bus, then mark it as completed (delete in the next cycle)
	{
		for(uint64_t j=0; j<R; j++)
		{
			if(s_queue[i].dest_reg == rb[j].dest_reg && s_queue[i].dest_tag == rb[j].tag && s_queue[i].filled_bit==1) 
			{
				s_queue[i].completed = 1;
				break;
			}
		}
	}

}

/**
 * Subroutine for cleaning up any outstanding instructions and calculating overall statistics
 * such as average IPC, average fire rate etc.
 * XXX: You're responsible for completing this routine
 *
 * @p_stats Pointer to the statistics structure
 */
void complete_proc(proc_stats_t *p_stats) 
{
	//Printing cycle by cylce behavior of each instruction in the trace
	cout<<"INST	"<<"FETCH	"<<"DISP	"<<"SCHED	"<<"EXEC	"<<"STATE	"<<endl;
	for(uint64_t i=1;i<=100000;i++)
	{
		cout<<i<<"	"<<out[i].fetched<<"	"<<out[i].fetched+1<<"	"<<out[i].dispatched+1<<"	"<<out[i].scheduled+1<<"	"<<out[i].status_update<<endl;
	}
	cout<<endl;
	p_stats->retired_instruction = number_of_retires;
	p_stats->cycle_count = cycle_count;
	p_stats->max_disp_size = maximum_dispatch_size;
	p_stats->avg_disp_size = total_disp_queue_size/(double)cycle_count;
	p_stats->avg_inst_retired = number_of_retires/(double)cycle_count;
	p_stats->avg_inst_fired = number_of_retires/(double)cycle_count;
	delete[] s_queue;
	delete[] f_unit;
	delete[] rb;
	delete[] reg_file;
	delete[] out;	
}
