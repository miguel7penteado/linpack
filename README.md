# linpack
Este repositório contém uma implementação do benchmark Linpack obtido de [netlib](http://www.netlib.org/benchmark/linpackc.new). Ele foi levemente modificado para não esperar nenhuma entrada de teclado.

Para executar o benchmark:
- `sudo docker run -it --rm elreyes/linpack`

Os resultados do benchmark serão impressos no console.

Aqui estão alguns resultados de exemplo do meu sistema com o tamanho de array padrão de 200:

```
$ sudo docker run -it --rm elreyes/linpack
Memory required:  315K.


LINPACK benchmark, Double precision.
Machine precision:  15 digits.
Array size 200 X 200.
Average rolled and unrolled performance:

    Reps Time(s) DGEFA   DGESL  OVERHEAD    KFLOPS
----------------------------------------------------
    2048   0.71  78.02%   2.89%  19.09%  4917805.960
    4096   1.43  78.02%   2.88%  19.11%  4870465.484
    8192   2.86  78.00%   2.87%  19.13%  4863107.385
   16384   5.70  77.98%   2.87%  19.14%  4882638.725
   32768  11.41  77.98%   2.87%  19.14%  4878150.976
```

É possível personalizar o tamanho do array usado pelo código de benchmark (o tamanho é 200 por padrão). Para fazer isso, defina uma variável de ambiente ao executar o contêiner:

`sudo docker run -it --rm -e LINPACK_ARRAY_SIZE=600 elreyes/linpack`

### Construindo a imagem localmente

Alternativamente, você pode clonar este repositório e construir sua própria imagem Docker localmente do zero. Para fazer isso, siga estas etapas:

- Clone este repositório (requer que o git esteja instalado em seu sistema... Google para seu SO/distro específico)
  - `git clone https://github.com/miguel7penteado/linpack.git`

- Construir a imagem do Docker
  - `cd linpack`
  - `sudo docker build -t linpack ./`
  
- Execute o benchmark
  - `sudo docker run -it --rm linpack`
  
referência:
[linpack](https://github.com/ereyes01/linpack.git)