`timescale 1ns / 1ps

module Newton_Raphson# (parameter iter = 11)
            (input clk, input reset, output done);

reg signed [31:0] x1 = 32'b00111111100000000000000000000000, x2 = 32'b01000000000000000000000000000000, x3 = 32'b01000000010000000000000000000000;
wire signed [31:0] J0,J1,J2,J3,J4,J5,J6,J7,J8;
wire signed [31:0] f1_1,f1_2,f1_3,f1_4,f1_5,f1_6,f1_int;
wire signed [31:0] f2_1,f2_2,f2_3,f2_4,f2_5,f2_6,f2_int;
wire signed [31:0] f3_1,f3_2,f3_3,f3_4,f3_5,f3_6,f3_7,f3_int;
wire signed [31:0] J4_1,J4_2,J8_1,J8_2,J8_3,J8_4;
wire signed [31:0] det_J;
wire signed [31:0] J00,J11,J22,J33,J44,J55,J66,J77,J88;
wire signed [31:0] J11_1,J11_2,J44_1,J44_2,J55_1,J55_2,J77_1,J77_2,J88_1,J88_2,J22_1,J22_2;
wire signed [31:0] det_1,det_2,det_3,det_4,det_5,det_6,det_7,det_8,det_9,det_10,det_11,det_12,det_13,det_14,det_15,det_16,det_13_int;
wire signed [31:0] f1,f2,f3;
wire signed [31:0] i0,i1,i2;
wire signed [31:0] i0_1,i0_2,i0_3,i0_4,i0_5,i1_1,i1_2,i1_3,i1_4,i1_5,i2_1,i2_2,i2_3,i2_4,i2_5;
wire signed [31:0] ans0,ans1,ans2;
reg [3:0] count;
///////////////////////////////////// f1 = x1*x1 - 2*x1 + x2*x2 - x3 + 1///////////////////////////////
Multiplication M1(x1,x1,,,,f1_1);
Addition_Subtraction A1(x1,x1,0,,f1_2);
Multiplication M2(x2,x2,,,,f1_3);
Addition_Subtraction A2(f1_1,f1_2,1,,f1_4);
Addition_Subtraction A3(f1_4,f1_3,0,,f1_5);
Addition_Subtraction A4(f1_5,x3,1,,f1_6);
Addition_Subtraction A5(f1_6,32'b00111111100000000000000000000000,0,,f1);

//////////////////////////////////////////f2 = x1*x2*x2 - x1 - 3*x2 + x2*x3 + 2 ////////////////////////////////

Multiplication M3(x1,f1_3,,,,f2_1);
Multiplication M4(32'b01000000010000000000000000000000,x2,,,,f2_2);
Multiplication M5(x2,x3,,,,f2_3);
Addition_Subtraction A6(f2_1,x1,1,,f2_4);
Addition_Subtraction A7(f2_4,f2_2,1,,f2_5);
Addition_Subtraction A8(f2_5,f2_3,0,,f2_6);
Addition_Subtraction A9(f2_6,32'b01000000000000000000000000000000,0,,f2);

/////////////////////////////////////// f3 = x1*x3*x3 - 3*x3 + x2*x3*x3 + x1*x2 ///////////////////////////////////

Multiplication M6(x3,x3,,,,f3_1);
Multiplication M7(x1,f3_1,,,,f3_2);
Multiplication M8(x3,32'b01000000010000000000000000000000,,,,f3_3);
Multiplication M9(f3_1,x2,,,,f3_4);
Multiplication M10(x1,x2,,,,f3_5);
Addition_Subtraction A10(f3_2,f3_3,1,,f3_6);
Addition_Subtraction A11(f3_6,f3_4,0,,f3_7);
Addition_Subtraction A12(f3_7,f3_5,0,,f3);


////////////////////////////////////////// Calculate Jacobian Matrix  //////////////////////////////////////////////////////////////

Addition_Subtraction A13(f1_2,32'b01000000000000000000000000000000,1,,J0);          // J0 = 2*x1-2;
Addition_Subtraction A14(x2,x2,0,,J1);                                              //assign J1 = 2*x2;
assign J2 = 32'b10111111100000000000000000000000;                                   //assign J2 = -1;
Addition_Subtraction A15(f1_3,32'b00111111100000000000000000000000,1,,J3);          //assign J3 = x2*x2-1;
Multiplication M11(32'b01000000000000000000000000000000,f3_5,,,,J4_1);
Addition_Subtraction A16(x3,J4_1,0,,J4_2);
Addition_Subtraction A17(J4_2,32'b01000000010000000000000000000000,1,,J4);          //assign J4 = x3+2*x1*x2-3;
assign J5 = x2;                                                                     //assign J5 = x2;
Addition_Subtraction A18(f3_1,J5,0,,J6);                                            //assign J6 = x3*x3+x2;
Addition_Subtraction A19(f3_1,x1,0,,J7);                                            //assign J7 = x3*x3+x1;
Multiplication M12(x1,x3,,,,J8_1);                      //  x1*x3
Addition_Subtraction A20(J8_1,J8_1,0,,J8_2);            //  2*x1*x3
Addition_Subtraction A21(f2_3,f2_3,0,,J8_3);            //  2*x2*x3
Addition_Subtraction A22(J8_2,J8_3,0,,J8_4);            //  2*x1*x3 + 2*x2*x3
Addition_Subtraction A23(J8_4,32'b01000000010000000000000000000000,1,,J8);          //assign J8 = 2*x1*x3+2*x2*x3-3;

////////////////////////////////////////Det of J /////////////////////////////////////////
//assign det_J = J0*J4*J8 - J0*J5*J7 - J1*J3*J8 + J1*J5*J6 + J2*J3*J7 - J2*J4*J6;

Multiplication M13(J4,J8,,,,det_1); 
Multiplication M14(J0,det_1,,,,det_2);
Multiplication M15(J5,J7,,,,det_3); 
Multiplication M16(J0,det_3,,,,det_4);
Multiplication M17(J3,J8,,,,det_5); 
Multiplication M18(J1,det_5,,,,det_6);
Multiplication M19(J5,J6,,,,det_7); 
Multiplication M20(J1,det_7,,,,det_8);
Multiplication M21(J3,J7,,,,det_9); 
Multiplication M22(J2,det_9,,,,det_10);
Multiplication M23(J4,J6,,,,det_11); 
Multiplication M24(J2,det_11,,,,det_12);

Addition_Subtraction A24(det_2,det_4,1,,det_13);
Addition_Subtraction A25(det_8,det_6,1,,det_14);
Addition_Subtraction A26(det_10,det_12,1,,det_15);
Addition_Subtraction A27(det_13,det_14,0,,det_16);
Addition_Subtraction A28(det_16,det_15,0,,det_J);

///////////////////////Inverse jacobian Matrix //////////////////////////

Addition_Subtraction A29(det_1,det_3,1,,J00);           //assign J00 = (J4*J8-J5*J7);
Multiplication M25(J1,J8,,,,J11_1); 
Multiplication M26(J2,J7,,,,J11_2); 
Addition_Subtraction A30(J11_2,J11_1,1,,J11);           //assign J11 = (-J1*J8+J2*J7);
Multiplication M27(J1,J5,,,,J22_1); 
Multiplication M28(J2,J4,,,,J22_2); 
Addition_Subtraction A31(J22_1,J22_2,1,,J22);           // assign J22 = (J1*J5-J4*J2);
Addition_Subtraction A32(det_7,det_5,1,,J33);           // assign J33 = (-J3*J8+J5*J6);
Multiplication M29(J0,J8,,,,J44_1); 
Multiplication M30(J6,J2,,,,J44_2); 
Addition_Subtraction A33(J44_1,J44_2,1,,J44);           // assign J44 = (J0*J8-J6*J2);
Multiplication M31(J0,J5,,,,J55_1); 
Multiplication M32(J3,J2,,,,J55_2); 
Addition_Subtraction A34(J55_2,J55_1,1,,J55);           // assign J55 = (-J0*J5+J2*J3);
Addition_Subtraction A35(det_9,det_11,1,,J66);           //assign J66 = J3*J7-J4*J6;
Multiplication M33(J0,J7,,,,J77_1); 
Multiplication M34(J1,J6,,,,J77_2); 
Addition_Subtraction A36(J77_2,J77_1,1,,J77);           // assign J77 = (-J0*J7+J1*J6);
Multiplication M35(J0,J4,,,,J88_1); 
Multiplication M36(J1,J3,,,,J88_2); 
Addition_Subtraction A37(J88_1,J88_2,1,,J88);            // assign J88 = (J4*J0-J1*J3);

Multiplication M37(f1,J00,,,,i0_1); 
Multiplication M38(f2,J11,,,,i0_2); 
Multiplication M39(f3,J22,,,,i0_3); 
Addition_Subtraction A38(i0_1,i0_2,0,,i0_4);
Addition_Subtraction A39(i0_3,i0_4,0,,i0_5);
Division d1(i0_5,det_J,,i0);                            // assign i0 = (f1*J00 + f2*J11 + f3*J22)/det_J;

Multiplication M40(f1,J33,,,,i1_1); 
Multiplication M41(f2,J44,,,,i1_2); 
Multiplication M42(f3,J55,,,,i1_3); 
Addition_Subtraction A40(i1_1,i1_2,0,,i1_4);
Addition_Subtraction A41(i1_3,i1_4,0,,i1_5);
Division d2(i1_5,det_J,,i1);                            // assign i1 = (f1*J33 + f2*J44 + f3*J55)/det_J;

Multiplication M43(f1,J66,,,,i2_1); 
Multiplication M44(f2,J77,,,,i2_2); 
Multiplication M45(f3,J88,,,,i2_3); 
Addition_Subtraction A42(i2_1,i2_2,0,,i2_4);
Addition_Subtraction A43(i2_3,i2_4,0,,i2_5);
Division d3(i2_5,det_J,,i2);                            // assign i2 = (f1*J66 + f2*J77 + f3*J88)/det_J;


Addition_Subtraction A44(x1,i0,1,,ans0);                //assign ans0 = x1 -  i0;
Addition_Subtraction A45(x2,i1,1,,ans1);                //assign ans1 = x2 -  i1;
Addition_Subtraction A46(x3,i2,1,,ans2);                //assign ans2 = x3 -  i2;

always @(posedge clk) begin
    if(count < iter & count != 0) begin
        x1 <= ans0;
        x2 <= ans1;
        x3 <= ans2;
    end
    else begin
        x1 <= x1;
        x2 <= x2;
        x3 <= x3;
    end
end
always @(posedge clk) begin
    if(reset) begin
        count <= 0;
    end
    else if (count<iter)
        count <= count + 1;
end

assign done = count==(iter-1) ? 1'b1 : 1'b0;

endmodule
