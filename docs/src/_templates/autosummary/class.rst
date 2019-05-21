{{ fullname | escape | underline}}

.. currentmodule:: {{ module }}

.. autoclass:: {{ objname }}

   {% if attributes %}
      .. rubric:: Attributes

      .. autosummary::
         :toctree:
      {% for item in attributes %}
         {{ name }}.{{ item }}
      {% endfor %}
   {% endif %}

   {% if methods %}
      .. rubric:: Methods

      .. autosummary::
         :toctree:
      {% for item in methods %}
         {{ name }}.{{ item }}
      {% endfor %}
   {% endif %}