`timescale 1ns / 1ps

module Newton_Raphson_TB( );

reg clk,reset;
parameter iter = 10;

Newton_Raphson#(.iter(iter))
     uut(.clk(clk),.reset(reset));

always #5 clk= ~clk;

initial begin
    $monitor("x1 =  %h , x2 = %h , x3 = %h", uut.x1,uut.x2,uut.x3);
    clk = 0;
    reset = 1;
    #10
    reset = 0;
end

endmodule
