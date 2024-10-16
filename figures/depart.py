from manim import *

config.background_color = WHITE

class Depart(Scene):
    def construct(self):
        labeli = MathTex(r"i").scale(0.5)
        labelim1 = MathTex(r"i-1").scale(0.5)
        labelip1 = MathTex(r"i+1").scale(0.5)
        labelj = MathTex(r"j=").scale(0.5)
        labelis = VGroup(labeli, labelip1.shift(UP), labelim1.shift(DOWN), labelj.shift(DOWN*2+RIGHT*0.25)).set_color(BLACK)

        I_1_square = Square(side_length=2, fill_color=BLUE, fill_opacity=1).scale(0.3)
        I_1_sqr_col = VGroup(I_1_square, I_1_square.copy().shift(UP), I_1_square.copy().shift(DOWN))
        labelkm1 = MathTex(r"k-1", color=BLUE).move_to(I_1_sqr_col).shift(DOWN*2).scale(0.5)

        I_1_dots = MathTex(r"\cdots", color=BLUE).move_to(I_1_sqr_col).shift(LEFT)
        I_1_sqr_col_1 = I_1_sqr_col.copy().move_to(I_1_dots).shift(LEFT)
        label1 = MathTex(r"1", color=BLUE).move_to(I_1_sqr_col_1).shift(DOWN*2).scale(0.5)
        labelI_1 = MathTex(r"I_1", color=BLUE).move_to(I_1_sqr_col).shift(UP*2).shift(LEFT)

        I_1 = VGroup(I_1_sqr_col, I_1_dots, I_1_sqr_col_1, label1, labelkm1, labelI_1)


        I_2_circle = Circle(radius=1, color=GREEN, fill_color=GREEN, fill_opacity=1).scale(0.3)
        I_2_circle_tri = VGroup(I_2_circle, I_2_circle.copy().shift(UP), I_2_circle.copy().shift(UR))
        labelk = MathTex(r"k", color=GREEN).move_to(I_2_circle).shift(DOWN*2).scale(0.5)
        labelI_2 = MathTex(r"I_2", color=GREEN).move_to(I_2_circle).shift(UP*2).shift(RIGHT*0.5)
        I_2 = VGroup(I_2_circle_tri, labelk, labelI_2)

        I_3_triangle = Triangle(color=RED, fill_color=RED, fill_opacity=1).shift(DOWN*0.25).scale(0.4)
        I_3_triangle_sk = VGroup(I_3_triangle, I_3_triangle.copy().shift(UR), I_3_triangle.copy().shift(DL))
        labelkp1 = MathTex(r"k+1", color=RED).move_to(I_3_triangle).shift(DOWN*2).scale(0.5)
        I_3_dots = MathTex(r"\cdots", color=RED).move_to(I_3_triangle_sk).shift(RIGHT)
        I_3_triangle_sk_2i_1 = I_3_triangle_sk.copy().move_to(I_3_dots).shift(RIGHT)
        label2i_1 = MathTex(r"m-1", color=RED).move_to(I_3_triangle_sk_2i_1).shift(DOWN*2).scale(0.5)
        labelI_3 = MathTex(r"I_3", color=RED).move_to(I_3_dots).shift(UP*2)

        I_3 = VGroup(I_3_triangle_sk, labelkp1, I_3_dots, I_3_triangle_sk_2i_1, label2i_1, labelI_3)

        I_4_circle = Circle(radius=1, color=GREEN_D, fill_color=GREEN_D, fill_opacity=1).scale(0.3)
        I_4_circle_tri = VGroup(I_4_circle, I_4_circle.copy().shift(DOWN), I_4_circle.copy().shift(DL))
        label2i = MathTex(r"m", color=GREEN_D).move_to(I_4_circle).shift(DOWN*2).scale(0.5)
        labelI_4 = MathTex(r"I_4", color=GREEN_D).move_to(I_4_circle).shift(UP*2).shift(LEFT*0.5)
        I_4 = VGroup(I_4_circle_tri, label2i, labelI_4)

        I_5_square = Square(side_length=2, fill_color=BLUE_D, fill_opacity=1).scale(0.3)
        I_5_sqr_col = VGroup(I_5_square, I_5_square.copy().shift(UP), I_5_square.copy().shift(DOWN))
        label2ip1 = MathTex(r"m+1", color=BLUE_D).move_to(I_5_sqr_col).shift(DOWN*2).scale(0.5)
        I_5_dots = MathTex(r"\cdots", color=BLUE_D).move_to(I_5_sqr_col).shift(RIGHT)
        I_5_sqr_col_2i_1 = I_5_sqr_col.copy().move_to(I_5_dots).shift(RIGHT)
        label2N = MathTex(r"2N", color=BLUE_D).move_to(I_5_sqr_col_2i_1).shift(DOWN*2).scale(0.5)
        labelI_5 = MathTex(r"I_5", color=BLUE_D).move_to(I_5_sqr_col).shift(UP*2).shift(RIGHT)
        I_5 = VGroup(I_5_sqr_col, label2ip1, labelI_5, I_5_dots, I_5_sqr_col_2i_1, label2N)

        self.add(I_1.shift(LEFT*3), I_2.shift(LEFT*2), I_3.shift(LEFT*1), labelis.shift(LEFT*6), I_4.shift(RIGHT*2), I_5.shift(RIGHT*3))