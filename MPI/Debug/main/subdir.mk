################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../main/bayeskit.c \
../main/bayeskitLV2.c \
../main/bayeskitMPI.c \
../main/bayeskitMPILV2.c \
../main/mpihw.c \
../main/mpisr.c \
../main/mpitest.c 

OBJS += \
./main/bayeskit.o \
./main/bayeskitLV2.o \
./main/bayeskitMPI.o \
./main/bayeskitMPILV2.o \
./main/mpihw.o \
./main/mpisr.o \
./main/mpitest.o 

C_DEPS += \
./main/bayeskit.d \
./main/bayeskitLV2.d \
./main/bayeskitMPI.d \
./main/bayeskitMPILV2.d \
./main/mpihw.d \
./main/mpisr.d \
./main/mpitest.d 


# Each subdirectory must supply rules for building sources it contributes
main/%.o: ../main/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -I/home/htihe9/code/mpi/sprng2.0/include -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


