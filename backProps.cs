using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Windows.Forms;

namespace BackForm
{
    public class backProp
    {
        public int R; /*層の数*/
        public int N1; /*入力層のユニット数*/
        public int NR; /*出力層のユニット数*/
        public int[] n; /*各層のユニット数*/
        public int MAX_N; /*層中の最大のユニット。LayerRecを作るために*/

        public double ETA_W;
        public double ETA_THETA;
        public int P;
        public double EPSILON;
        public double LOOP;
        public double[,] T;
        public int columns;
        public int rows;
        public double[] pattern;

        public double[,] tempo;
        public int[] extra_lessons;
        public int EL;
        public ConvertNR nr;
        public double[,] result;
        public double[,] maxs_mins;

        //それぞれの層ごとに、ユニットのデータを作る
        public struct LayerRec
        {

            public double[,] w;
            public double[] theta;
            public double[] output;
            public double[] delta;
            public void Reset(int MAX_N) {
                w = new double[MAX_N, MAX_N];
                theta = new double[MAX_N];
                output = new double[MAX_N];
                delta = new double[MAX_N];
            }

        }

        public LayerRec[] L;


        /*
         * コンストラクタ
         * 全ての層の全てのユニットの結合加重と閾値に0-1までの少数をセットする。
         */
        public backProp(double p_loop, double p_eta, double p_epsilon, int hidden_L, int hidden_U,int output_u,int EL, double[,] list)
        {




            columns = list.GetLength(0); // 行
            rows = list.GetLength(1); // 列 つまり入力ユニット+出力ユニット

            //教師信号のパターンの数
            P = columns;
            pattern = new double[P];

            //許容誤差。誤差関数がこれに近づくようにする
            EPSILON = p_epsilon;

            R = 2+ hidden_L ; /*層の数*/
            N1 = rows - output_u; /*入力層のユニット数*/
            NR = rows - N1; /*出力層のユニット数*/



            n = new int[R];
            n[0] = N1;
            int q;
            for (q = 1; q < R-1; q++) {
                n[q] = hidden_U;
            }
            n[q] = NR; /*各層のユニット数*/

            
            MAX_N = rows; // 層中の最大のユニット。LayerRecを作るために。適当にrows

            //学習定数
            ETA_W = p_eta;
            ETA_THETA = p_eta;

            //学習回数
            LOOP = p_loop;

            result =new double[P,NR];

            //チェック用
            tempo = new double[P, P];

            //正規化
            T = new double[columns, rows];
            this.nr = new ConvertNR(ref list);
            T = nr.Normalization();
            this.maxs_mins = nr.maxs_mins;


            //補習回数記録配列
            extra_lessons = new int[P];

            //補習回数
            this.EL = EL;
            

            L = new LayerRec[R];



            int i, j, k;
            Random cR = new Random();

            //↓これは入力層に入力<->出力するだけの層、つまり L[0].output[i] の配列だけ確保している
            //入力層に先立つ第0層扱い
            L[0].Reset(MAX_N);

            //結合加重をセット。[-1,1]の範囲で乱数を発生させてる
            for (k=1;k<R;k++){
                L[k].Reset(MAX_N);
                for (i = 0; i < n[k]; i++)
                {
                    for (j = 0; j < n[k - 1]; j++)
                        L[k].w[i, j] = 2*cR.NextDouble()-1.0;
                    L[k].theta[i] = 2 * cR.NextDouble() - 1.0;
                }
            }
        }


        //入力→出力
        public void Forward (int p){
            int i, j, k;
            double u;

            k = 0;
            for (i = 0; i < n[k]; i++)
                L[k].output[i] = T[p, i];
            for(k=1;k<R;k++)
                for (i = 0; i < n[k]; i++)
                {
                    u = 0.0;
                    for (j = 0; j < n[k - 1]; j++)
                        u += L[k].w[i, j] * L[k - 1].output[j];
                    u -= L[k].theta[i];
                    L[k].output[i] = sigmo(u);
                }
        }

        //出力→入力
        public void Backward(int p) {
            int i, k, l;
            double delta;

            k = R - 1;//一番最後の層=出力層
            //↓出力層の全てのユニットに対し
            for (i = 0; i < n[k]; i++)
                L[k].delta[i] = -(T[p, i+N1] - L[k].output[i]) * f_dash(L[k].output[i]);
            //中間層の計算
            for (k = R - 2; k >= 1; k--)
                for (i = 0; i < n[k]; i++)
                {
                    delta = 0.0;
                    for (l = 0; l < n[k + 1]; l++)
                        delta += L[k + 1].w[l, i] * L[k + 1].delta[l];
                    L[k].delta[i] = delta * f_dash(L[k].output[i]);
                }
                
        }
        //学習
        public void Learn()
        {
            int i, j, k;
            for(k=R-1;k>=1;k--)
                for (i = 0; i < n[k]; i++)
                {
                    for (j = 0; j < n[k - 1]; j++)
                        L[k].w[i, j] += -ETA_W * L[k].delta[i] * L[k - 1].output[j];
                    L[k].theta[i] += ETA_THETA * L[k].delta[i];
                }
        }


        //誤差を測る + 補習機能
        //優秀な生徒は卒業
        public int ExtraLessons(int p) {
            int i,k, flag;
            k = R - 1;
            double exE = new double();

            flag = 0;

            for (i = 0; i < n[k]; i++)
            {
                exE += Math.Abs(T[p, i + N1] - L[k].output[i]);
            }
            pattern[p] = exE;
            if (exE < EPSILON) {
                extra_lessons[p]++;
                flag = 1;
                if (EL <= extra_lessons[p]) {
                    flag = 2;
                }
            }

            return flag;
        }


        // 全て学習
        // 戻り値はループが何回できたか回数
        public double LearnAll() {
            int p, loop,clear,flag;



            flag = 0;
            loop = 0;
            for(;;){
                clear = 0;
                for (p = 0; p < P; p++) {
                    Forward(p);
                    if(EL != 0)
                        flag = ExtraLessons(p);
                    if (flag >= 1) {
                        clear++;
                    }
                    if(!(flag >= 2))
                    {
                        Backward(p);
                        Learn();
                    }

                }
                if (clear >= columns)
                {
                    break;
                }

                if (loop > LOOP)
                {
                    break;
                }

                loop++;
            
            }
            Test1();
            return loop-1;
        }


        //予測値の抽出1
        public void Test1() {
            int p, i, j, k, l;
            double u;

            for (p = 0; p < P; p++)
            {

                k = 0;
                for (i = 0; i < n[k]; i++)
                    L[k].output[i] = T[p, i];
                for (k = 1; k < R; k++)
                    for (i = 0; i < n[k]; i++)
                    {
                        u = 0.0;
                        for (j = 0; j < n[k - 1]; j++)
                            u += L[k].w[i, j] * L[k - 1].output[j];
                        u -= L[k].theta[i];
                        L[k].output[i] = sigmo(u);
                    }
                for (l = 0; l < n[R-1]; l++)
                    result[p, l] = L[R-1].output[l];
            }
            this.nr.Specialization(ref result);
        
        }


        //感度分析用
        public double[,] AnalyzeSense() {
            int i,k,l,j,h;
            double[] based = new double[NR];
            double[,] sense = new double[N1, NR];
            double u,difference;
            

            //全て最小値のときの出力(base)
            k = 0;
            for (i = 0; i < n[k]; i++)
                L[k].output[i] = 0;
            for (k = 1; k < R; k++)
                for (i = 0; i < n[k]; i++)
                {
                    u = 0.0;
                    for (j = 0; j < n[k - 1]; j++)
                        u += L[k].w[i, j] * L[k - 1].output[j];
                    u -= L[k].theta[i];
                    L[k].output[i] = sigmo(u);
                }
            for (l = 0; l < n[R - 1]; l++)
                based[l] = L[R - 1].output[l];



            //それぞれの最大値での計算
            for (j = 0; j < N1; j++)
            {
                
                k = 0;
                for (i = 0; i < N1; i++)
                {
                    if (j == i)
                    {
                        difference = maxs_mins[j, 0] - maxs_mins[j, 1];
                        if (difference == 0)
                            L[k].output[i] = 0.0;
                        else
                            L[k].output[i] = 1.0;

                    }
                    else
                        L[k].output[i] = 0.0;
                }
                for (k = 1; k < R; k++)
                {
                    for (i = 0; i < n[k]; i++)
                    {
                        u = 0.0;
                        for (h = 0; h < n[k - 1]; h++)
                            u += L[k].w[i, h] * L[k - 1].output[h];
                        u -= L[k].theta[i];

                        L[k].output[i] = sigmo(u);
                    }
                }
                for (l = 0; l < n[R - 1]; l++)
                    sense[j, l] = L[R - 1].output[l] - based[l];
             }
            return sense;
        
        
        }


        //テスト 調べたいpの値を入れて調べる
        //予測用。
        public double[,] Test(ref double[,] list)
        {
            int p,i, j, k,l;
            double u,difference;
            double[,] tmp;
            double[,] tmp2;

            int columns = list.GetLength(0); // 行
            int rows = list.GetLength(1); // 列 つまり入力ユニット+出力ユニット

            tmp = new double[columns, rows];


            //最大値・最小値から正規化する。
            for (i = 0; i < rows; i++)
            {
                difference = maxs_mins[i, 0] - maxs_mins[i, 1];
                for (l = 0; l < columns; l++)
                {
                    tmp[l, i] = (list[l, i] - maxs_mins[i, 1]) / difference;
                    if (difference == 0) { tmp[l, i] = 0; }
                }
            }




            double[,] answer = new double[columns, NR];
            

            //フォワード処理
            for (p = 0; p < columns; p++)
            {

                k = 0;
                for (i = 0; i < n[k]; i++)
                    L[k].output[i] = tmp[p, i];
                for (k = 1; k < R; k++)
                    for (i = 0; i < n[k]; i++)
                    {
                        u = 0.0;
                        for (j = 0; j < n[k - 1]; j++)
                            u += L[k].w[i, j] * L[k - 1].output[j];
                        u -= L[k].theta[i];
                        L[k].output[i] = sigmo(u);
                    }
                for (l = 0; l < n[R - 1]; l++)
                    answer[p, l] = L[R - 1].output[l];
            }

            tmp2 = new double[columns, NR];

            //特殊化した数字を代入。入力ユニットの数の分はいらないのでN1を足して
            //出力ユニットの分の最大値・最小値で計算
            for (i = 0; i < NR; i++)
            {
                difference = maxs_mins[i + rows, 0] - maxs_mins[i + rows, 1];
                for (l = 0; l < columns; l++)
                {
                    tmp2[l, i] = (answer[l, i] * difference) + maxs_mins[i + rows, 1];
                }

            }


            return tmp2;
        }





        //シグモイド関数を計算して返すだけ
        public double sigmo(double a)
        {
            return 1.0 / (1.0 + Math.Exp(-(a)));
        }

        //シグモイドを微分した形
        public double f_dash(double a)
        {
            return (a * (1.0 - a));
        }

    }
    public class ConvertNR {
        public double[,] maxs_mins;
        public double[,] list;
        public int columns;
        public int rows;

        public ConvertNR(ref double[,] list) {
            int i, l;
            this.list = list;


            columns = list.GetLength(0); // 行
            rows = list.GetLength(1); // 列 つまりこの場合入力ユニット

            this.maxs_mins = new double[rows,2];

            double[] tmp = new double[columns];


            //その列の最大値・最小値を
            for (i = 0; i < rows; i++)
            {
                for (l = 0; l < columns; l++)
                {
                    tmp[l] = list[l, i];
                }
                maxs_mins[i, 0] = tmp.Max();
                maxs_mins[i, 1] = tmp.Min();
            }
        
        }
        //与えられたリストから正規化
        //変換後リストを戻り値
        public double[,] Normalization()
        {
            int i, l;
            double difference;
            double[,] tmp;

            tmp = new double[columns, rows];


            //最大値・最小値から正規化する。
            for (i = 0; i < rows; i++)
            {
                difference = maxs_mins[i, 0] - maxs_mins[i, 1];
                for (l = 0; l < columns; l++)
                {
                    tmp[l, i] = (list[l, i] - maxs_mins[i, 1]) / difference;
                    if (difference == 0) { tmp[l, i] = 0; }
                }
            }
            list = tmp;
            return tmp;

        }
        //与えられたリストと最大値・最小値リストから元の値に戻す
        //変換後リストはリファレンスで変更
        public void Specialization(ref double[,] exanswer)
        {
            int i, l;
            double difference;
            double[,] tmp;



            int columns = exanswer.GetLength(0); // アンサー配列の行
            int rows = exanswer.GetLength(1); // アンサー配列の列 

            tmp = new double[columns, rows];

            //特殊化した数字を代入。入力ユニットの数の分はいらないのでN1を足して
            //出力ユニットの分の最大値・最小値で計算
            for (i = 0; i < rows; i++)
            {
                difference = maxs_mins[i + this.rows - rows, 0] - maxs_mins[i + this.rows - rows, 1];
                for (l = 0; l < columns; l++)
                {
                    tmp[l, i] = (exanswer[l, i] * difference) + maxs_mins[i + this.rows - rows, 1];
                }

            }
            exanswer = tmp;


        }



    }

}
